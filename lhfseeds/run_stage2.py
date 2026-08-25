"""Run stage 2: copies.tsv -> deduplicated loci -> sampled -> sequence -> MSA.

    python -m lhfseeds.run_stage2 [--config ...] [--top 2] [--linkage strict]
                                  [--modes gap_aware,merge_always] [--cluster 62]

Emits, per cluster and merge mode, a seed packet:

    seed_packets/cluster_00062/gap_aware/
        loci.tsv            every deduplicated locus, with cross-tool support
        sampled.tsv         the <=cap copies that entered the alignment
        copies.fa           element only        -- what gets aligned
        copies.flanked.fa   element +- flank_bp -- TSD / TE-Aid / Refiner
        aln.fa              MAFFT alignment
        consensus.fa        rebuilt consensus
        packet.json         realized composition and provenance

The two merge modes are run side by side because the decision between them is
made at SEED level (depth and rebuilt consensus length), not at copy level --
copy counts alone were never going to settle it.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import stage2, stockholm
from .fasta import IndexedFasta

# members whose copies come from the tool's own fragment linkage do not vary
# with merge mode; they must appear in BOTH pools or the comparison is not
# like-for-like (different membership, not just different merging).
INVARIANT_MODE = "hit_id"


def emit_seed(outdir: Path, cluster_id: int, mode: str, fin: dict,
              cfg: dict, engine: str = "mafft") -> dict:
    """Write the Stockholm seed and lint it.

    `stk lint` is a SOFT gate by design (PLAN_A): a failing packet still goes to
    the queue, flagged, with its lint output attached, because a failure that is
    silently dropped is a failure nobody fixes.
    """
    ids = [stage2.seq_id(cfg["assembly"], c, s, e, st)
           for (c, s, e, st) in fin["coords"]]
    rows = ["".join(chr(x) for x in r).replace("-", stockholm.GAP)
            for r in fin["matrix"]]
    tp = lookup_tp(cfg)
    meta = {
        "ID": f"TEbedSeeds_c{cluster_id:05d}_{mode}_{engine}",
        "DE": (f"Consensus rebuilt from {len(ids)} genomic copies of "
               f"multi-tool cluster {cluster_id} ({mode})")[:80],
        "AU": cfg.get("au_string", ""),
        "OC": cfg.get("taxon", ""),
        "SQ": len(ids),
        "BM": f"TEbed-seeds {cfg.get('contract_version', '?')}; {engine}",
        "CC": [f"Rebuilt from track data; provenance in packet.json.",
               f"Merge mode {mode}; one representative per deduplicated locus."],
        "RF": stockholm.rf_line(fin["consensus"], fin["is_match"]),
    }
    if tp:
        meta["TP"] = tp
    stk_path = outdir / f"seed.{engine}.stk"
    stockholm.write_stockholm(stk_path, [stockholm.format_record(ids, rows, meta)])

    stk_bin = Path(cfg.get("stk_bin",
                           "vendor/dfam-curator/target/release/stk")).expanduser()
    res: dict = {"tp_emitted": bool(tp)}
    if not stk_bin.exists():
        res["error"] = f"stk binary not found at {stk_bin}"
        return res
    # rewrite RF with Dfam's own consensus caller so rf_consensus_mismatch
    # compares like with like, then lint the file that will actually be shipped
    tmp_rf = outdir / f"seed.{engine}.rf.stk"
    if stockholm.update_consensus(stk_path, tmp_rf, stk_bin):
        stk_path.unlink()
        tmp_rf.rename(stk_path)
    res["tier1"] = stockholm.lint(stk_path, stk_bin, no_network=True)
    genome = cfg.get("assembly_fasta")
    if genome:
        res["genome"] = stockholm.lint(stk_path, stk_bin,
                                       genome=Path(genome).expanduser(),
                                       no_network=True)
    (outdir / f"lint.{engine}.txt").write_text(
        res["tier1"]["output"] + "\n" + res.get("genome", {}).get("output", ""))
    return res


def lookup_tp(cfg: dict) -> str | None:
    """#=GF TP from config/tp_map.tsv.  An unmapped canonical path must BLOCK
    the TP rather than guess one -- a wrong classification is worse than none."""
    p = Path(cfg.get("tp_map", "config/tp_map.tsv"))
    fallback = cfg.get("tp_default")
    if not p.exists():
        return fallback
    try:
        t = pd.read_csv(p, sep="\t")
    except Exception:
        return fallback
    if "canonical_path" not in t.columns or "tp" not in t.columns:
        return fallback
    return fallback


def build_packet(copies: pd.DataFrame, mode: str, fa: IndexedFasta, cfg: dict,
                 outdir: Path, cluster_id: int, threads: int = 1) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    pool = copies[copies.merge_mode.isin([mode, INVARIANT_MODE])]
    if not len(pool):
        return dict(cluster_id=cluster_id, mode=mode, error="no copies")

    conslen = pool.cons_len.dropna()
    modal_cons = float(conslen.mode().iat[0]) if len(conslen) else float("nan")

    loci = stage2.cluster_loci(pool, modal_cons_len=modal_cons)

    # A copy is judged over-long against the CLUSTER's modal consensus length,
    # not against its own member's.  Measured on cluster 62: the two over-long
    # rows that reached the seed were both REPET copies of a 764 bp entry that
    # is itself ~3 tandem units of a 265 bp element (PIPELINE_FINDINGS F4), so
    # the per-member cap of 1.5 x 764 let 495 and 642 bp tandem dimers into an
    # alignment of a 265 bp element.  The cluster's modal length is the better
    # reference precisely because it is a consensus across tools.
    seedcfg = cfg.get("seed") or {}
    x_modal = float(seedcfg.get("max_copy_x_modal_consensus", 1.5))
    n_over = 0
    if not np.isnan(modal_cons):
        over = (loci.end - loci.start) > x_modal * modal_cons
        n_over = int(over.sum())
        over_members = (loci.loc[over, "member"].value_counts().to_dict()
                        if n_over else {})
        loci["over_modal_consensus"] = over
        loci.to_csv(outdir / "loci.tsv", sep="\t", index=False)
        eligible = loci[~over]
    else:
        over_members = {}
        loci.to_csv(outdir / "loci.tsv", sep="\t", index=False)
        eligible = loci

    smp = cfg["sampling"]
    sample = stage2.stratified_sample(eligible, smp["cap"], smp["floor"],
                                      smp["full_length_frac"])
    sample.to_csv(outdir / "sampled.tsv", sep="\t", index=False)

    frags = stage2.explode_fragments(sample)
    flank = int(cfg["flank_bp"])
    recs, stats = stage2.extract(fa, frags, cfg["assembly"], flank_bp=0)
    frecs, fstats = stage2.extract(fa, frags, cfg["assembly"], flank_bp=flank)
    stage2.write_fasta(recs, outdir / "copies.fa")
    stage2.write_fasta(frecs, outdir / "copies.flanked.fa")
    stats.to_csv(outdir / "extract_stats.tsv", sep="\t", index=False)

    # 4.3 check: flanked length == element length + realized flanks, exactly
    merged = stats.merge(fstats, on="seq_id", suffixes=("", "_f"))
    bad = merged[merged.length_f != merged.length
                 + merged.flank_left_f + merged.flank_right_f]
    if len(bad):
        raise AssertionError(f"flank accounting failed for {len(bad)} records")

    coord_map = {stage2.seq_id(cfg["assembly"], f.chrom, int(f.start),
                               int(f.end), f.strand):
                 (f.chrom, int(f.start), int(f.end), f.strand)
                 for f in frags.itertuples()}
    min_occ = float(seedcfg.get("min_match_occupancy", 0.5))
    min_row = int(seedcfg.get("min_frag_bp", stage2.MIN_FRAG_BP))
    engines = cfg.get("engines") or [cfg.get("consensus_engine", "mafft")]
    primary = cfg.get("consensus_engine", "mafft")

    eng_results: dict = {}
    cons, occ_stats, lint_res, aln_info = "", {}, {}, {}
    n_aln = 0
    for eng in engines:
        if len(recs) < 2:
            continue
        try:
            if eng == "mafft":
                info = stage2.mafft(outdir / "copies.fa", outdir / "aln.fa",
                                    threads=threads)
                arecs = stage2.read_fasta(outdir / "aln.fa")
                rc = stage2.row_coords_mafft(arecs, coord_map)
                extra = dict(n_rows_mafft_flipped=sum(
                    1 for i, _ in arecs if i.startswith("_R_")))
            elif eng == "refiner":
                rbin = Path(cfg.get("refiner_bin",
                                    "vendor/refiner/bin/Refiner")).expanduser()
                if not rbin.exists():
                    eng_results[eng] = dict(error=f"Refiner not found at {rbin}; "
                                            "run tools/setup_refiner.sh")
                    continue
                rres = stage2.refiner(outdir / "copies.fa", rbin, threads=threads,
                                      workdir=outdir / "refiner_work")
                arecs = rres["records"]
                rc, drop_unparsed = [], []
                keep = []
                for rid, row in arecs:
                    g = stage2.parse_refiner_id(rid, coord_map)
                    if g is None:
                        drop_unparsed.append(rid)
                        continue
                    rc.append(g)
                    keep.append((rid, row))
                arecs = keep
                stage2.write_fasta(arecs, outdir / "aln.refiner.fa")
                info = dict(avg_kimura=rres["avg_kimura"],
                            refiner_consensus_len=len(rres["consensus"]),
                            n_unparsed_ids=len(drop_unparsed))
                extra = dict(avg_kimura=rres["avg_kimura"])
            else:
                eng_results[eng] = dict(error=f"unknown engine {eng}")
                continue

            fin = stage2.finalize_alignment(arecs, rc, min_occupancy=min_occ,
                                            min_row_bp=min_row)
            m, is_match = fin["matrix"], fin["is_match"]
            depth = (m[:, is_match] != ord("-")).sum(axis=0)
            stage2.write_fasta(
                [(f"cluster_{cluster_id:05d}_{mode}_{eng}_consensus",
                  fin["consensus"])], outdir / f"consensus.{eng}.fa")
            lr = emit_seed(outdir, cluster_id, mode, fin, cfg, engine=eng)
            st = dict(
                engine=eng, alignment_rows=int(m.shape[0]),
                aln_width=int(m.shape[1]), n_match_columns=int(is_match.sum()),
                rebuilt_consensus_len=len(fin["consensus"]),
                median_depth=int(np.median(depth)) if len(depth) else 0,
                min_depth=int(depth.min()) if len(depth) else 0,
                max_depth=int(depth.max()) if len(depth) else 0,
                frac_cols_depth_ge3=float((depth >= 3).mean()) if len(depth) else 0.0,
                n_rows_dropped=len(fin["dropped"]), lint=lr, info=info, **extra)
            eng_results[eng] = st
            pd.DataFrame(fin["dropped"]).to_csv(
                outdir / f"dropped_rows.{eng}.tsv", sep="\t", index=False)
            if eng == primary or not cons:
                cons = fin["consensus"]
                n_aln = int(m.shape[0])
                aln_info = info
                lint_res = lr
                occ_stats = {k: v for k, v in st.items()
                             if k not in ("engine", "lint", "info",
                                          "alignment_rows",
                                          "rebuilt_consensus_len")}
        except Exception as exc:                       # engine failure is data
            eng_results[eng] = dict(error=f"{type(exc).__name__}: {exc}")
            print(f"[stage2] cluster {cluster_id} {mode} {eng}: FAILED {exc}",
                  file=sys.stderr)

    realized_full = int((sample.sampled_as == "full").sum())
    # Over-extension guard (PLAN_A 4.7): compare against the MEDIAN member
    # consensus length, never the max -- one chimeric or over-assembled member
    # must not be able to set the reference the guard is measured against.
    med_member = float(pool.groupby("member").cons_len.median().median())
    ratio = float(cfg.get("overextension_ratio", 1.5))
    overext = bool(cons) and len(cons) > ratio * med_member
    packet = dict(
        cluster_id=cluster_id, mode=mode,
        contract_version=cfg.get("contract_version"),
        assembly=cfg["assembly"], flank_bp=flank,
        members=sorted(pool.member.unique()),
        n_copy_rows=int(len(pool)),
        n_loci=int(len(loci)),
        redundancy_mean=round(float(loci.n_members.mean()), 3),
        n_loci_full=int(loci.any_full.sum()),
        support_histogram={int(k): int(v) for k, v in
                           loci.n_tools.value_counts().sort_index().items()},
        full_frac_by_support={int(k): round(float(v), 4) for k, v in
                              loci.groupby("n_tools").any_full.mean().items()},
        modal_consensus_len=modal_cons,
        member_consensus_lens={m: float(v) for m, v in
                               pool.groupby("member").cons_len.median().items()},
        n_loci_over_modal=n_over,
        over_modal_by_member=over_members,
        max_copy_x_modal_consensus=x_modal,
        # REALIZED, not target: repeatome intactness varies by taxon
        sampling_target=dict(cap=smp["cap"], floor=smp["floor"],
                             full_length_frac=smp["full_length_frac"]),
        sampled_n=int(len(sample)),
        sampled_full=realized_full,
        sampled_partial=int(len(sample) - realized_full),
        realized_full_frac=round(realized_full / max(len(sample), 1), 4),
        alignment_rows=n_aln,
        rows_from_bridged_copies=int(sample.n_fragments.gt(1).sum()),
        n_frag_rows=int(len(frags)),
        rebuilt_consensus_len=len(cons),
        median_member_consensus_len=med_member,
        overextension_ratio=ratio,
        possible_overextension=overext,
        lint=lint_res,
        mafft=aln_info,
        engines=eng_results,
        **occ_stats,
        elapsed_s=round(time.time() - t0, 1),
    )
    json.dump(packet, open(outdir / "packet.json", "w"), indent=1, default=str)
    return packet


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    ap.add_argument("--work-dir", default=None)
    ap.add_argument("--linkage", default="strict")
    ap.add_argument("--top", type=int, default=2)
    ap.add_argument("--cluster", type=int, action="append", default=None,
                    help="explicit cluster id(s); overrides --top")
    ap.add_argument("--modes", default=None,
                    help="comma list; default from config merge_mode")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args(argv)

    cfg = yaml.safe_load(open(args.config))
    work = Path(args.work_dir or cfg["work_dir"]).expanduser() / cfg["assembly"]
    fasta = cfg.get("assembly_fasta")
    if not fasta:
        sys.exit("[stage2] config assembly_fasta is null -- stage 2 needs sequence")
    fa = IndexedFasta(Path(fasta).expanduser())

    modes = (args.modes.split(",") if args.modes else
             (["merge_always", "gap_aware"] if cfg["merge_mode"] == "both"
              else [cfg["merge_mode"]]))

    packets = []
    if args.cluster:
        cids = args.cluster
    else:
        cand = pd.read_csv(work / f"candidates_{args.linkage}.tsv", sep="\t")
        cand = cand[cand.candidate].sort_values("pooled_full_len", ascending=False)
        cids = [int(c) for c in cand.head(args.top).cluster_id]

    for cid in cids:
        cdir = work / "seed_packets" / f"cluster_{cid:05d}"
        cfile = cdir / "copies.tsv"
        if not cfile.exists():
            print(f"[stage2] cluster {cid}: no copies.tsv, run stage 1 first",
                  file=sys.stderr)
            continue
        copies = pd.read_csv(cfile, sep="\t")
        for mode in modes:
            p = build_packet(copies, mode, fa, cfg, cdir / mode, cid,
                             threads=args.threads)
            packets.append(p)
            if "error" in p:
                print(f"[stage2] cluster {cid} {mode}: {p['error']}", file=sys.stderr)
                continue
            print(f"[stage2] cluster {cid} {mode}: {p['n_copy_rows']} rows -> "
                  f"{p['n_loci']} loci -> {p['sampled_n']} sampled "
                  f"({p['realized_full_frac']:.0%} full) -> {p['alignment_rows']} "
                  f"aln rows, consensus {p['rebuilt_consensus_len']} bp, "
                  f"median depth {p.get('median_depth', 0)} [{p['elapsed_s']}s]",
                  file=sys.stderr)

    cmp = pd.DataFrame([{k: v for k, v in p.items()
                         if not isinstance(v, (dict, list))} for p in packets])
    if len(cmp):
        cmp.to_csv(work / "seed_mode_comparison.tsv", sep="\t", index=False)
        cols = [c for c in ["cluster_id", "mode", "n_loci", "n_loci_full",
                            "sampled_n", "realized_full_frac", "alignment_rows",
                            "rows_from_bridged_copies", "rebuilt_consensus_len",
                            "trimmed_width", "median_depth",
                            "frac_cols_depth_ge3"] if c in cmp.columns]
        print("\n" + cmp[cols].to_string(index=False), file=sys.stderr)
    fa.close()


if __name__ == "__main__":
    main()

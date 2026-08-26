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
              cfg: dict, engine: str = "mafft",
              tp_info: tuple = (None, "", "unmapped", None),
              fa: IndexedFasta | None = None,
              cluster_key: str | None = None,
              kimura: float | None = None,
              model: dict | None = None) -> dict:
    """Write the Stockholm seed and lint it.

    `stk lint` is a SOFT gate by design (PLAN_A): a failing packet still goes to
    the queue, flagged, with its lint output attached, because a failure that is
    silently dropped is a failure nobody fixes.
    """
    ids = [stage2.seq_id(cfg["assembly"], c, s, e, st)
           for (c, s, e, st) in fin["coords"]]
    rows = ["".join(chr(x) for x in r).replace("-", stockholm.GAP)
            for r in fin["matrix"]]
    tp, tp_ver, tp_status, cpath = tp_info
    # A two-model cluster emits two seeds; without the model in the ID they
    # collide and `stk lint` raises duplicate_id when the batch is concatenated.
    mtag = ""
    if model and model.get("n_models", 1) > 1:
        mtag = f"_m{model['index']}x{int(round(model['centre']))}"
    meta = {
        "ID": f"TEbedSeeds_c{cluster_id:05d}_{mode}{mtag}_{engine}",
        "DE": (f"Consensus rebuilt from {len(ids)} genomic copies of "
               f"multi-tool cluster {cluster_id} ({mode}"
               + (f"; {int(round(model['centre']))} bp length mode" if mtag else "")
               + ")")[:80],
        "AU": cfg.get("au_string", ""),
        "OC": cfg.get("taxon", ""),
        "SQ": len(ids),
        "BM": f"TEbed-seeds {cfg.get('contract_version', '?')}; {engine}",
        # SE is capped at 80 characters (se_too_long, ERROR)
        "SE": (f"TEbed-seeds {mode}; cluster {cluster_key or cluster_id}; "
               f"{cfg['assembly']}")[:80],
        "CC": [f"Rebuilt from track data; provenance in packet.json.",
               f"Merge mode {mode}; one representative per deduplicated locus.",]
               + ([f"Length model {model['index'] + 1} of {model['n_models']} "
                   f"for this cluster, centre {int(round(model['centre']))} bp, "
                   f"supported by {model['n_tools']} tools "
                   f"({','.join(model['tools'])}). The other model(s) describe "
                   f"the same family at a different length -- for an LTR family "
                   f"the solo LTR and the full element."] if mtag else [])
               + [
               f"Classification path {cpath or 'unresolved'}"
               + (f"; TP scheme {tp_ver}" if tp_ver else "")],
        "RF": stockholm.rf_line(fin["consensus"], fin["is_match"]),
    }
    if tp:
        meta["TP"] = tp
    if kimura is not None and kimura == kimura:      # not NaN
        meta["KD"] = f"{kimura:.2f}"
    stk_path = outdir / f"seed.{engine}.stk"
    stockholm.write_stockholm(stk_path, [stockholm.format_record(ids, rows, meta)])

    stk_bin = Path(cfg.get("stk_bin",
                           "vendor/dfam-curator/target/release/stk")).expanduser()
    res: dict = {"tp_emitted": bool(tp), "tp": tp,
                 "tp_scheme_version": tp_ver, "tp_status": tp_status,
                 "canonical_path": cpath}
    if not stk_bin.exists():
        res["error"] = f"stk binary not found at {stk_bin}"
        return res
    # rewrite RF with Dfam's own consensus caller so rf_consensus_mismatch
    # compares like with like, then lint the file that will actually be shipped
    tmp_rf = outdir / f"seed.{engine}.rf.stk"
    if stockholm.update_consensus(stk_path, tmp_rf, stk_bin):
        stk_path.unlink()
        tmp_rf.rename(stk_path)
    # Our own coordinate check, because lint's is advisory: --genome reports
    # every coordinate problem as WARN and still exits 0, and an unparseable
    # identifier is dropped from validation silently.
    if fa is not None:
        res["coord_check"] = stage2.verify_rows_against_genome(
            fin["coords"], rows, fa)
    res["tier1"] = stockholm.lint(stk_path, stk_bin, no_network=True)
    genome = cfg.get("assembly_fasta")
    if genome:
        res["genome"] = stockholm.lint(stk_path, stk_bin,
                                       genome=Path(genome).expanduser(),
                                       no_network=True)
        coord_codes = {"seq_coord_invalid", "seq_coord_fixed",
                       "seq_id_not_in_ref"}
        res["genome_coord_warnings"] = sum(
            n for c, n in res["genome"]["codes"].items()
            if c.split(":", 1)[-1] in coord_codes)
    (outdir / f"lint.{engine}.txt").write_text(
        res["tier1"]["output"] + "\n" + res.get("genome", {}).get("output", ""))
    return res


def load_tp_map(cfg: dict) -> dict:
    """canonical_path -> (tp, scheme_version, status)."""
    p = Path(cfg.get("tp_map", "config/tp_map.tsv"))
    if not p.exists():
        return {}
    t = pd.read_csv(p, sep="\t", comment="#")
    if "canonical_path" not in t.columns or "tp" not in t.columns:
        return {}
    return {r.canonical_path: (r.tp, getattr(r, "scheme_version", ""),
                               getattr(r, "status", ""))
            for r in t.itertuples()}


def lookup_tp(tp_map: dict, canonical_path, cfg: dict) -> tuple:
    """#=GF TP for a cluster's canonical path.

    An unmapped path must BLOCK the TP rather than guess one: `tp_unknown` is a
    lint ERROR, and a seed carrying a confidently WRONG classification is worse
    than one carrying none -- the first is silently absorbed into Dfam, the
    second is visibly incomplete. The packet records which scheme version was
    used, because the Dfam scheme is being reconciled with Repbase and every
    seed must say what it was classified against.
    """
    if canonical_path and canonical_path in tp_map:
        tp, ver, status = tp_map[canonical_path]
        return tp, ver, status, canonical_path
    return None, "", "unmapped", canonical_path


def build_packet(copies: pd.DataFrame, mode: str, fa: IndexedFasta, cfg: dict,
                 outdir: Path, cluster_id: int, threads: int = 1,
                 tp_info: tuple = (None, "", "unmapped", None),
                 cluster_key: str | None = None) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    pool = copies[copies.merge_mode.isin([mode, INVARIANT_MODE])]
    if not len(pool):
        return dict(cluster_id=cluster_id, mode=mode, error="no copies")

    seedcfg = cfg.get("seed") or {}
    x_modal = float(seedcfg.get("max_copy_x_modal_consensus", 1.5))
    split_ratio = float(seedcfg.get("mode_split_ratio", 2.0))
    min_mode_tools = int(seedcfg.get("min_tools_per_mode", 2))

    # A first pass with no cap, only to get the deduplicated loci. The modal
    # length MUST be taken over loci, not over pool copy rows: a member with
    # many redundant annotations otherwise decides the cluster's length. On
    # cluster 1183 the pool-row mode is 423 bp and the locus mode is 7400 bp,
    # and the 423 answer deleted every full-length locus.
    loci = stage2.cluster_loci(pool)
    # Only MEASURED lengths vote on the modes. A member whose cons_len_source is
    # `clustermate` or `deconvolved` is repeating a cluster-mate's opinion, not
    # supplying independent evidence, so counting its tool inflates the mode's
    # apparent cross-tool support -- and a #=GF CC line claiming N tools support
    # a length would then ship a provenance claim that is not true. Measured on
    # cluster 1183: edta contributes 423 and 7400 bp, but EDTA's own library
    # entries are 244 and 272 bp and its BED carries NA consensus coordinates.
    meas = loci[loci.cons_len_source == "bed16"] if "cons_len_source" in loci else loci
    member_lens = (meas.groupby("member").cons_len.median().dropna().to_dict()
                   if "cons_len" in meas.columns else {})
    if not member_lens and "cons_len" in loci.columns:   # nothing measured
        member_lens = loci.groupby("member").cons_len.median().dropna().to_dict()
    modes = stage2.length_modes(member_lens, split_ratio, min_mode_tools)
    credible = [m for m in modes if m["credible"]]
    if not credible and modes:
        credible = [max(modes, key=lambda m: m["n_members"])]
    centres = [m["centre"] for m in credible]
    modal_cons = float(centres[0]) if len(centres) == 1 else (
        float(np.median(centres)) if centres else float("nan"))

    # Each locus is capped against the mode IT belongs to. A 7.4 kb locus in a
    # cluster whose credible modes are 423 bp and 7367 bp is a full element,
    # not a 17-fold tandem multimer of the solo LTR, and capping it against
    # 423 was what threw away the evidence.
    spans = (loci.end - loci.start).to_numpy()
    if centres:
        midx = stage2.assign_mode(spans, centres)
        ref = np.array(centres)[midx]
        over = spans > x_modal * ref
        loci["mode_index"] = midx
        loci["mode_centre"] = ref
    else:
        over = np.zeros(len(loci), dtype=bool)
        loci["mode_index"] = 0
        loci["mode_centre"] = np.nan
    loci["over_modal_consensus"] = over
    n_over = int(over.sum())
    over_members = (loci.loc[over, "member"].value_counts().to_dict()
                    if n_over else {})
    loci.to_csv(outdir / "loci.tsv", sep="\t", index=False)
    eligible = loci[~over]

    # One model per credible length mode. With a single mode this is exactly
    # the previous behaviour and writes flat into outdir; with two it writes
    # model0_.../model1_... subdirectories and the cluster-level packet.json
    # lists them, so a curator sees the solo LTR and the full element as the
    # pair they are rather than one silently winning.
    models = []
    for i, m in enumerate(credible):
        models.append(dict(index=i, centre=m["centre"], n_tools=m["n_tools"],
                           tools=m["tools"], members=m["members"],
                           n_models=len(credible)))
    if len(models) <= 1:
        only = models[0] if models else None
        return _build_model(eligible, loci, pool, mode, fa, cfg, outdir,
                            cluster_id, threads, tp_info, cluster_key,
                            seedcfg, x_modal, modal_cons, n_over,
                            over_members, modes, only, t0)

    built = []
    for m in models:
        sub = eligible[eligible.mode_index == m["index"]]
        sdir = outdir / f"model{m['index']}_{int(round(m['centre']))}bp"
        if len(sub) < 2:
            built.append(dict(model_index=m["index"], centre=m["centre"],
                              error=f"only {len(sub)} loci in this mode"))
            continue
        built.append(_build_model(sub, loci, pool, mode, fa, cfg, sdir,
                                  cluster_id, threads, tp_info, cluster_key,
                                  seedcfg, x_modal, modal_cons, n_over,
                                  over_members, modes, m, t0))
    packet = dict(
        cluster_id=cluster_id, cluster_key=cluster_key, mode=mode,
        multi_model=True, n_models=len(models),
        cluster_modal_consensus_len=modal_cons,
        length_modes=[{k: v for k, v in mm.items() if k != "members"}
                      for mm in modes],
        models=[{k: v for k, v in b.items()
                 if k in ("model_index", "model_centre", "model_n_tools",
                          "model_tools", "n_loci", "sampled_n",
                          "rebuilt_consensus_len", "error")} for b in built],
        model_dirs=[f"model{m['index']}_{int(round(m['centre']))}bp"
                    for m in models],
        n_loci=int(len(loci)), n_loci_over_modal=n_over,
        elapsed_s=round(time.time() - t0, 1),
    )
    json.dump(packet, open(outdir / "packet.json", "w"), indent=1, default=str)
    return packet


def _build_model(eligible: pd.DataFrame, loci: pd.DataFrame,
                 pool: pd.DataFrame, mode: str, fa: IndexedFasta,
                 cfg: dict, outdir: Path, cluster_id: int,
                 threads: int, tp_info: tuple, cluster_key,
                 seedcfg: dict, x_modal: float, modal_cons: float,
                 n_over: int, over_members: dict, modes: list,
                 model: dict | None, t0: float) -> dict:
    """Build one consensus model from an eligible set of loci.

    A cluster with two credible length modes (an LTR family with a solo
    LTR and a full element) yields two models, each built from its own
    loci. Splitting here rather than picking one length is the point of
    the change: picking one threw the other away.
    """
    outdir.mkdir(parents=True, exist_ok=True)
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
            # A seed with no match column has no consensus. It is structurally
            # valid Stockholm and lints CLEAN -- 184 such seeds shipped from the
            # first 200-cluster batch before this check existed -- so it must be
            # caught here or not at all.
            degenerate = int(is_match.sum()) == 0
            depth = (m[:, is_match] != ord("-")).sum(axis=0)
            stage2.write_fasta(
                [(f"cluster_{cluster_id:05d}_{mode}_{eng}_consensus",
                  fin["consensus"])], outdir / f"consensus.{eng}.fa")
            lr = emit_seed(outdir, cluster_id, mode, fin, cfg, engine=eng,
                           tp_info=tp_info, fa=fa,
                           cluster_key=cluster_key,
                           kimura=extra.get("avg_kimura"), model=model)
            st = dict(
                engine=eng, alignment_rows=int(m.shape[0]),
                aln_width=int(m.shape[1]), n_match_columns=int(is_match.sum()),
                rebuilt_consensus_len=len(fin["consensus"]),
                median_depth=int(np.median(depth)) if len(depth) else 0,
                min_depth=int(depth.min()) if len(depth) else 0,
                max_depth=int(depth.max()) if len(depth) else 0,
                frac_cols_depth_ge3=float((depth >= 3).mean()) if len(depth) else 0.0,
                n_rows_dropped=len(fin["dropped"]), degenerate=degenerate,
                lint=lr, info=info, **extra)
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
        cluster_id=cluster_id, cluster_key=cluster_key, mode=mode,
        contract_version=cfg.get("contract_version"),
        assembly=cfg["assembly"], flank_bp=flank,
        members=sorted(pool.member.unique()),
        n_copy_rows=int(len(pool)),
        n_loci=int(len(eligible)),
        n_loci_cluster=int(len(loci)),
        redundancy_mean=round(float(loci.n_members.mean()), 3),
        n_loci_full=int(loci.any_full.sum()),
        support_histogram={int(k): int(v) for k, v in
                           loci.n_tools.value_counts().sort_index().items()},
        full_frac_by_support={int(k): round(float(v), 4) for k, v in
                              loci.groupby("n_tools").any_full.mean().items()},
        canonical_path=tp_info[3], tp=tp_info[0],
        tp_scheme_version=tp_info[1], tp_status=tp_info[2],
        modal_consensus_len=(model["centre"] if model else modal_cons),
        cluster_modal_consensus_len=modal_cons,
        model_index=(model["index"] if model else 0),
        model_centre=(model["centre"] if model else modal_cons),
        model_n_tools=(model["n_tools"] if model else None),
        model_tools=(model["tools"] if model else None),
        model_members=(model["members"] if model else None),
        n_models=(model["n_models"] if model else 1),
        length_modes=[{k: v for k, v in m.items() if k != "members"}
                      for m in modes],
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
    ap.add_argument("--clusters-file", default=None,
                    help="TSV with a cluster_id column (tools/select_batch.py)")
    ap.add_argument("--modes", default=None,
                    help="comma list; default from config merge_mode")
    ap.add_argument("--threads", type=int, default=4,
                    help="threads handed to MAFFT within one packet")
    ap.add_argument("--workers", type=int, default=1,
                    help="packets built in parallel, as separate processes")
    ap.add_argument("--shard", default=None, metavar="I/N",
                    help="internal: process only clusters where index %% N == I")
    ap.add_argument("--skip-existing", action="store_true",
                    help="skip clusters whose packet.json already exists "
                         "for every requested mode")
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
    tp_map = load_tp_map(cfg)
    cand = pd.read_csv(work / f"candidates_{args.linkage}.tsv", sep="\t")
    path_by_cluster = (dict(zip(cand.cluster_id, cand.majority_path))
                       if "majority_path" in cand.columns else {})
    key_by_cluster = (dict(zip(cand.cluster_id, cand.cluster_key))
                      if "cluster_key" in cand.columns else {})
    if not path_by_cluster:
        print("[stage2] candidates table has no majority_path column -- re-run "
              "stage 0 to emit TP classifications", file=sys.stderr)
    if args.clusters_file:
        cids = [int(c) for c in
                pd.read_csv(args.clusters_file, sep="\t").cluster_id]
    elif args.cluster:
        cids = args.cluster
    else:
        sel = cand[cand.candidate].sort_values("pooled_full_len", ascending=False)
        cids = [int(c) for c in sel.head(args.top).cluster_id]

    if args.skip_existing:
        before = len(cids)
        cids = [c for c in cids
                if not all((work / "seed_packets" / f"cluster_{c:05d}" / m
                            / "packet.json").exists() for m in modes)]
        print(f"[stage2] --skip-existing: {before - len(cids)} already built",
              file=sys.stderr)

    # Parallel by process over disjoint packet directories, same as stage 1.
    if args.workers > 1 and args.shard is None:
        import subprocess
        base = [sys.executable, "-m", "lhfseeds.run_stage2",
                "--config", args.config, "--linkage", args.linkage,
                "--threads", str(max(1, args.threads // args.workers))]
        if args.work_dir:
            base += ["--work-dir", args.work_dir]
        if args.modes:
            base += ["--modes", args.modes]
        for c in cids:
            base += ["--cluster", str(int(c))]
        if args.skip_existing:
            base += ["--skip-existing"]
        procs = [subprocess.Popen(base + ["--shard", f"{i}/{args.workers}"])
                 for i in range(args.workers)]
        codes = [p.wait() for p in procs]
        bad = [i for i, c in enumerate(codes) if c != 0]
        if bad:
            sys.exit(f"[stage2] shard(s) {bad} failed with {codes}")
        print(f"[stage2] {args.workers} shards complete", file=sys.stderr)
        return

    if args.shard:
        i, n = (int(x) for x in args.shard.split("/"))
        cids = list(cids)[i::n]

    for cid in cids:
        cdir = work / "seed_packets" / f"cluster_{cid:05d}"
        cfile = cdir / "copies.tsv"
        if not cfile.exists():
            print(f"[stage2] cluster {cid}: no copies.tsv, run stage 1 first",
                  file=sys.stderr)
            continue
        copies = pd.read_csv(cfile, sep="\t")
        for mode in modes:
            cpath = path_by_cluster.get(cid)
            cpath = None if pd.isna(cpath) else cpath
            tp_info = lookup_tp(tp_map, cpath, cfg)
            p = build_packet(copies, mode, fa, cfg, cdir / mode, cid,
                             threads=args.threads, tp_info=tp_info,
                             cluster_key=key_by_cluster.get(cid))
            packets.append(p)
            if "error" in p:
                print(f"[stage2] cluster {cid} {mode}: {p['error']}", file=sys.stderr)
                continue
            if p.get("multi_model"):
                bits = " | ".join(
                    f"model{m['model_index']} ~{int(m['model_centre'])}bp -> "
                    f"{m.get('rebuilt_consensus_len', '?')} bp"
                    for m in p.get("models", []))
                print(f"[stage2] cluster {cid} {mode}: {p['n_loci']} loci, "
                      f"{p['n_models']} length models: {bits}", file=sys.stderr)
                continue
            print(f"[stage2] cluster {cid} {mode}: {p['n_copy_rows']} rows -> "
                  f"{p['n_loci']} loci -> {p['sampled_n']} sampled "
                  f"({p['realized_full_frac']:.0%} full) -> {p['alignment_rows']} "
                  f"aln rows, consensus {p['rebuilt_consensus_len']} bp, "
                  f"median depth {p.get('median_depth', 0)} [{p['elapsed_s']}s]",
                  file=sys.stderr)

    # Batch lint triage (PLAN_A 4.5). `stk lint` is a SOFT gate: a failing
    # packet is still queued, flagged, with its lint output attached. The point
    # of the triage table is that the failure CODES are counted across the
    # batch, so a systematic format error shows up as one row to fix rather
    # than as N packets quietly dropped.
    tri = []
    for p in packets:
        for eng, est in (p.get("engines") or {}).items():
            if "error" in est:
                tri.append(dict(cluster_id=p["cluster_id"], mode=p["mode"],
                                engine=eng, tier="engine", severity="ERROR",
                                code="engine_failed", n=1,
                                detail=est["error"][:200]))
                continue
            if est.get("degenerate"):
                tri.append(dict(cluster_id=p["cluster_id"], mode=p["mode"],
                                engine=eng, tier="build", severity="ERROR",
                                code="degenerate_alignment", n=1,
                                detail="no match column; consensus is empty"))
            for tier in ("tier1", "genome"):
                res = (est.get("lint") or {}).get(tier) or {}
                for code, n in (res.get("codes") or {}).items():
                    sev, _, name = code.partition(":")
                    tri.append(dict(cluster_id=p["cluster_id"], mode=p["mode"],
                                    engine=eng, tier=tier, severity=sev,
                                    code=name, n=n, detail=""))
    tri_df = pd.DataFrame(tri)
    if len(tri_df):
        sfx = f".shard{args.shard.replace('/', '_')}" if args.shard else ""
        tri_df.to_csv(work / f"lint_triage{sfx}.tsv", sep="\t", index=False)
        roll = (tri_df.groupby(["severity", "code"])
                      .agg(packets=("cluster_id", "size"), total=("n", "sum"))
                      .sort_values(["severity", "total"], ascending=[True, False]))
        print("\n[stage2] lint triage across the batch:", file=sys.stderr)
        print(roll.to_string(), file=sys.stderr)
        n_err = int(tri_df[tri_df.severity == "ERROR"].cluster_id.nunique())
        print(f"[stage2] packets with >=1 lint ERROR: {n_err} / {len(packets)} "
              f"(soft gate -- all are still queued, each with lint.<engine>.txt)",
              file=sys.stderr)

    cmp = pd.DataFrame([{k: v for k, v in p.items()
                         if not isinstance(v, (dict, list))} for p in packets])
    if len(cmp):
        sfx2 = f".shard{args.shard.replace('/', '_')}" if args.shard else ""
        cmp.to_csv(work / f"seed_mode_comparison{sfx2}.tsv", sep="\t", index=False)
        cols = [c for c in ["cluster_id", "mode", "n_loci", "n_loci_full",
                            "sampled_n", "realized_full_frac", "alignment_rows",
                            "rows_from_bridged_copies", "rebuilt_consensus_len",
                            "trimmed_width", "median_depth",
                            "frac_cols_depth_ge3"] if c in cmp.columns]
        print("\n" + cmp[cols].to_string(index=False), file=sys.stderr)
    fa.close()


if __name__ == "__main__":
    main()

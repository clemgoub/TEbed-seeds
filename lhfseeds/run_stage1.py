"""Run stage 1 on the top-N candidate clusters.

    python -m lhfseeds.run_stage1 [--config ...] [--top 5] [--linkage strict]

Coordinate-only prototype: emits copies.tsv per cluster (Smitten ids,
per-member boundary claims, both merge modes side by side) and a
merge-mode comparison table. Sequence extraction waits on assembly_fasta.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import stage1


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    ap.add_argument("--work-dir", default=None)
    ap.add_argument("--linkage", default="strict")
    ap.add_argument("--top", type=int, default=5)
    args = ap.parse_args(argv)

    cfg = yaml.safe_load(open(args.config))
    repo = Path(cfg["vgp_tebed_repo"]).expanduser()
    work = Path(args.work_dir or cfg["work_dir"]).expanduser() / cfg["assembly"]
    cand = pd.read_csv(work / f"candidates_{args.linkage}.tsv", sep="\t")
    cand = cand[cand.candidate].sort_values("pooled_full_len", ascending=False)
    print(f"[stage1] {int(cand.shape[0])} candidates; processing top {args.top}",
          file=sys.stderr)

    sz = pd.read_csv(repo / "data" / f"{cfg['assembly']}.chrom.sizes", sep="\t",
                     names=["chrom", "size"])
    chrom_sizes = dict(zip(sz.chrom, sz["size"]))
    CLUSTER_TOOLS = ["rm2", "edta", "pantera", "fastltr", "repet"]

    comp_rows = []
    for _, crow in cand.head(args.top).iterrows():
        cid = int(crow.cluster_id)
        members = crow.members.split(";")
        cdir = work / "seed_packets" / f"cluster_{cid:05d}"
        cdir.mkdir(parents=True, exist_ok=True)
        cluster_fams = {tuple(m.split(":", 1)) for m in members}

        # one paint per cluster, reused for every gap check (was: per-gap rescan)
        t_probe = time.time()
        _, foreign_cov = stage1.build_foreign_probe(
            repo, CLUSTER_TOOLS, chrom_sizes, cluster_fams)
        print(f"[stage1] cluster {cid}: foreign probe painted in "
              f"{time.time()-t_probe:.0f}s", file=sys.stderr)

        # pass 1: own consensus lengths (modal), so mates can be inherited from
        hits_by_member, conslen = {}, {}
        for m in members:
            tool, fam = m.split(":", 1)
            h = stage1.member_hits(repo, tool, fam)
            if not len(h):
                continue
            hits_by_member[m] = h
            conslen[m] = stage1.consensus_length(h)
        # --- over-assembly deconvolution (PIPELINE_FINDINGS F4) ---------------
        # A member's library consensus can be several tandem units of the real
        # element collapsed into one entry: cluster 62's REPET entries are
        # 764/765 bp where rm2, pantera and edta all say 264-269. Judging those
        # copies against their own consensus makes real copies look 71%
        # truncated -- measured, 0.74% near-full-length against a cross-tool
        # ground truth of 49.5%, an under-count of ~67x.
        #
        # The discriminator is the CROSS-TOOL ratio, not the shape of the
        # coordinates. A coordinate-only detector (length commensurability plus
        # block structure) was built and calibrated against a half-integer-k
        # decoy null, which no tandem over-assembly can produce: real 70 vs
        # decoy 74.7, ratio 0.94, p = 0.72 -- no measurable specificity. So the
        # rule here uses only what the cluster provides, which is exactly what
        # the cluster is for.
        #
        # Inherit the mate length; do NOT divide by a detected k. k is not
        # identified by coordinates, and of 18 flagged families with
        # length-carrying mates only 4 had L/mate ~ k -- for the other 14 the
        # consensus was already SHORTER than their mates', so dividing would
        # have compounded the error up to 5x.
        # The rule abstains unless the MATES AGREE WITH EACH OTHER. Inheriting
        # from cluster-mates can otherwise inherit an over-assembly: a cluster
        # with only two length-carrying members, one of them collapsed, has a
        # median that is already wrong. Requiring >=2 other members whose
        # lengths agree within `mate_spread_max` makes the reference a genuine
        # cross-tool consensus rather than a single opinion.
        deconvolved, abstained = {}, {}
        ratio_min = float(cfg.get("overassembly_min_ratio", 1.8))
        spread_max = float(cfg.get("overassembly_mate_spread_max", 1.3))
        lens = {m: v for m, v in conslen.items() if not np.isnan(v)}
        for m, own in lens.items():
            others = [v for k, v in lens.items() if k != m]
            if len(others) < 2:
                continue
            mate = float(np.median(others))
            spread = max(others) / max(min(others), 1e-9)
            if mate > 0 and own / mate >= ratio_min and spread > spread_max:
                abstained[m] = dict(own=own, mate=mate,
                                    mate_spread=round(spread, 2),
                                    reason="cluster-mates disagree with each "
                                           "other; cannot tell which is the unit")
                print(f"[stage1] cluster {cid}: {m} is {own/mate:.2f}x its "
                      f"mates but the mates spread {spread:.2f}x -- ABSTAINING",
                      file=sys.stderr)
                continue
            if mate > 0 and own / mate >= ratio_min:
                deconvolved[m] = dict(own=own, mate=mate,
                                      ratio=round(own / mate, 3),
                                      k_implied=round(own / mate, 1))
                conslen[m] = mate
                print(f"[stage1] cluster {cid}: {m} consensus {own:.0f} bp is "
                      f"{own/mate:.2f}x its cluster-mates ({mate:.0f} bp) -- "
                      f"treating as over-assembled, using the mate length",
                      file=sys.stderr)

        mate_lookup = stage1.build_clustermate_conslen(
            repo, cluster_fams, chrom_sizes, conslen)

        all_copies = []
        for m, hits in hits_by_member.items():
            tool = m.split(":", 1)[0]
            own_len = conslen[m]   # deconvolved above if over-assembled
            use_hitid = tool == "repet" and hits.hit_id.notna().any()
            for mode in (["hit_id"] if use_hitid else ["merge_always", "gap_aware"]):
                cp = stage1.merge_copies(
                    hits, mode, foreign_cov if mode == "gap_aware" else None,
                    cons_len=own_len)
                cp["member"] = m
                cp["merge_mode"] = mode
                if np.isnan(own_len):
                    # inherit per LOCUS from whichever cluster-mate covers it
                    cp["cons_len"] = [mate_lookup(r.chrom, r.start, r.end)
                                      for r in cp.itertuples()]
                    cp["cons_len_source"] = np.where(cp.cons_len.notna(),
                                                     "clustermate", "none")
                    # the runaway-chain cap could not run inside merge_copies
                    # (no length was known then); re-split now that it is
                    over = (cp.end - cp.start) > (
                        stage1.MAX_COPY_X_CONSENSUS * cp.cons_len)
                    if over.any():
                        kept, split = cp[~over.fillna(True)], cp[over.fillna(False)]
                        pieces = [kept]
                        for r in split.itertuples():
                            fs = [int(x) for x in str(r.frag_starts).split(";")]
                            fe = [int(x) for x in str(r.frag_ends).split(";")]
                            cap = stage1.MAX_COPY_X_CONSENSUS * r.cons_len
                            cur_s, cur_e, sub = [fs[0]], [fe[0]], []
                            for s, e in zip(fs[1:], fe[1:]):
                                if e - cur_s[0] > cap:
                                    sub.append((cur_s[:], cur_e[:]))
                                    cur_s, cur_e = [s], [e]
                                else:
                                    cur_s.append(s); cur_e.append(e)
                            sub.append((cur_s, cur_e))
                            for ss, ee in sub:
                                pieces.append(pd.DataFrame([dict(
                                    chrom=r.chrom, start=min(ss), end=max(ee),
                                    strand=r.strand, n_fragments=len(ss),
                                    frag_starts=";".join(map(str, ss)),
                                    frag_ends=";".join(map(str, ee)),
                                    div=r.div, cons_start=r.cons_start,
                                    cons_end=r.cons_end, member=m,
                                    merge_mode=mode, cons_len=r.cons_len,
                                    cons_len_source=r.cons_len_source)]))
                        cp = pd.concat(pieces, ignore_index=True)
                else:
                    cp["cons_len"] = own_len
                    cp["cons_len_source"] = ("deconvolved" if m in deconvolved
                                             else "bed16")
                # full-length needs per-copy consensus coords, which an
                # inheriting member does not have -- flag, do not fake
                cp["is_full"] = (
                    stage1.near_full_length(cp, own_len,
                                            span_only=m in deconvolved)
                    if not np.isnan(own_len)
                    else pd.Series(False, index=cp.index))
                cp["smitten_id"] = [stage1.smitten(cfg["assembly"], r)
                                    for r in cp.itertuples()]
                all_copies.append(cp)

        copies = pd.concat(all_copies, ignore_index=True)
        copies.to_csv(cdir / "copies.tsv", sep="\t", index=False)
        if deconvolved or abstained:
            pd.DataFrame(
                [dict(member=m, action="deconvolved", **v)
                 for m, v in deconvolved.items()]
                + [dict(member=m, action="abstained", **v)
                   for m, v in abstained.items()]
            ).to_csv(cdir / "deconvolved.tsv", sep="\t", index=False)
        for mode in ["merge_always", "gap_aware", "hit_id"]:
            sub = copies[copies.merge_mode == mode]
            if not len(sub):
                continue
            ln = sub.end - sub.start
            cl = sub.cons_len.replace(0, np.nan)
            comp_rows.append(dict(
                cluster_id=cid, mode=mode, n_copies=len(sub),
                n_full=int(sub.is_full.sum()),
                med_len=int(ln.median()), max_len=int(ln.max()),
                max_x_consensus=round(float((ln / cl).max()), 1),
                frac_over_1p5x=round(float(((ln / cl) > 1.5).mean()), 3),
                multi_frag=int((sub.n_fragments > 1).sum()),
            ))
        print(f"[stage1] cluster {cid}: {len(copies)} copy rows -> {cdir}", file=sys.stderr)

    comp = pd.DataFrame(comp_rows)
    comp.to_csv(work / "merge_mode_comparison.tsv", sep="\t", index=False)
    print(comp.to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()

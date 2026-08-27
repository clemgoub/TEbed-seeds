"""Build the LHF curator queue: size + protein-domain gate, then rank.

    PYTHONPATH=. .venv/bin/python tools/lhf_queue.py \
        --work-dir runs/GCA_951799975.1 --out runs/GCA_951799975.1/lhf

Replaces ranking by `pooled_full_len`, which is a COPY-NUMBER ranking and so
promotes exactly what low-hanging fruit should exclude: abundant short
non-autonomous elements. The first top-20 had a median rebuilt consensus of
469 bp with 19 of 20 under 1 kb.

Two gates, then a rank:

  size    per-order minimum consensus length (SINEs exempt but capped)
  domain  at least one CORE replication domain at >=50% HMM coverage,
          found with BATH against TE-Aid's curated TE-protein Pfam profiles
  rank    within tier A, by distinct core domains then near-full-length copies

Both halves are needed. Domain presence rises steeply with size (0% below
300 bp, 84.5% above 2.5 kb) so size alone is a decent proxy -- but short
domain-positive models exist and would be wrongly discarded, and oversized
domain-free models exist and would be wrongly promoted; those are the chimera
and rDNA suspects a size-only filter puts at the top of the queue.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lhfseeds import domains as dom          # noqa: E402

BATH = Path("/Users/clementgoubert/Documents/BATH/opt/bin/bathsearch")
PFAM = Path("/Users/clementgoubert/Documents/TE-Aid/dev-data/protein-cache/"
            "07532cf00f976128/pfam.bhmm")


def collect_models(work: Path, engine: str) -> pd.DataFrame:
    """One row per built model, single- and multi-model layouts alike."""
    rows = []
    pats = ["seed_packets/cluster_*/*/packet.json",
            "seed_packets/cluster_*/*/model*/packet.json"]
    for pat in pats:
        for f in sorted(work.glob(pat)):
            try:
                p = json.load(open(f))
            except Exception:
                continue
            if "error" in p or not p.get("engines"):
                continue
            e = (p["engines"].get(engine) or {})
            if "error" in e or not e.get("rebuilt_consensus_len"):
                continue
            cons = f.parent / f"consensus.{engine}.fa"
            if not cons.exists():
                continue
            model = f.parent.name if f.parent.name.startswith("model") else "single"
            sid = f"c{p['cluster_id']:05d}_{p['mode']}_{model}"
            rows.append(dict(
                sid=sid, cluster=p["cluster_id"], cluster_key=p.get("cluster_key"),
                mode=p["mode"], model=model, fasta=str(cons),
                bp=e["rebuilt_consensus_len"], aln_rows=e.get("alignment_rows", 0),
                median_depth=e.get("median_depth", 0),
                n_full=p.get("n_loci_full", 0), n_loci=p.get("n_loci", 0),
                n_tools=len(p.get("members") or []) and
                        len({m.split(":")[0] for m in p["members"]}),
                tp=p.get("tp"), canonical_path=p.get("canonical_path"),
            ))
    d = pd.DataFrame(rows)
    if len(d):
        d["order"] = [dom.order_of(r.canonical_path, r.tp) for r in d.itertuples()]
    return d


def write_all_consensi(models: pd.DataFrame, out_fa: Path) -> None:
    with open(out_fa, "w") as fh:
        for r in models.itertuples():
            seq = "".join(l.strip() for l in open(r.fasta)
                          if not l.startswith(">"))
            if not seq:
                continue
            fh.write(f">{r.sid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", default="runs/GCA_951799975.1")
    ap.add_argument("--out", default=None)
    ap.add_argument("--engine", default="refiner")
    ap.add_argument("--bath", default=str(BATH))
    ap.add_argument("--profiles", default=str(PFAM))
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--cpu", type=int, default=8)
    args = ap.parse_args(argv)

    work = Path(args.work_dir)
    out = Path(args.out or work / "lhf")
    out.mkdir(parents=True, exist_ok=True)

    models = collect_models(work, args.engine)
    if not len(models):
        raise SystemExit(f"no {args.engine} models under {work}/seed_packets")
    print(f"models: {len(models)}  clusters: {models.cluster.nunique()}  "
          f"engine: {args.engine}")

    fa = out / "all_consensi.fa"
    write_all_consensi(models, fa)
    total_bp = int(models.bp.sum())
    print(f"scanning {total_bp/1e6:.2f} Mb with BATH ...")
    info = dom.run_bath(fa, out / "domains.tbl", Path(args.bath),
                        Path(args.profiles), cpu=args.cpu)
    if info["returncode"] != 0:
        print(f"  BATH failed: {info['stderr_tail'][:400]}")
    hits = dom.parse_bath_tbl(out / "domains.tbl")
    print(f"  hits: {len(hits)}  sequences with a hit: {hits.seq.nunique() if len(hits) else 0}")

    summ = dom.summarise_domains(hits)
    d = models.merge(summ, left_on="sid", right_on="seq", how="left")
    for c in ("n_dom", "n_core"):
        d[c] = d[c].fillna(0).astype(int)
    d["best_cov"] = d.best_cov.fillna(0.0)
    d["core_domains"] = d.core_domains.fillna("")
    d["domains"] = d.domains.fillna("")
    d["min_len"] = d.order.map(dom.MIN_LEN)
    d["tier"] = [dom.assign_tier(r) for _, r in d.iterrows()]

    d.sort_values(["tier", "n_core", "bp"], ascending=[True, False, False]) \
     .to_csv(out / "lhf_tiers.tsv", sep="\t", index=False)

    print("\n=== tiers ===")
    print(d.tier.value_counts().sort_index().to_string())
    print("\n=== tier by order ===")
    print(pd.crosstab(d.order, d.tier.str[0]).to_string())

    a = dom.rank_tier_a(d)
    a.to_csv(out / "lhf_tierA_ranked.tsv", sep="\t", index=False)
    print(f"\n=== tier A, top {args.top} (ranked by core domains, then full-length copies) ===")
    cols = ["sid", "order", "tp", "bp", "n_core", "core_domains", "best_cov",
            "n_full", "aln_rows", "median_depth"]
    show = a.head(args.top)[cols].copy()
    show["core_domains"] = show.core_domains.str.slice(0, 46)
    show["best_cov"] = show.best_cov.round(2)
    print(show.to_string(index=False) if len(show) else "  (empty)")
    print(f"\n  tier A total: {len(a)}   written to {out}/lhf_tierA_ranked.tsv")


if __name__ == "__main__":
    main()

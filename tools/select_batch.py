"""Choose a stratified batch of candidate clusters to build.

    PYTHONPATH=. .venv/bin/python tools/select_batch.py --n 200 --out batch.tsv

Why not just take the top of `pooled_full_len`: the 12 clusters built so far
are all in the **top 2%** by copy number and 9 of 12 are 4-tool. That is
precisely the stratum where tools agree, so it is the worst possible sample for
testing whether engine disagreement tracks classification depth (F8b).

The sample is stratified on two axes that matter for that test:

  * `path_depth` — how far the member-weighted canonical classification path
    resolves. 3 = stops at Class I/II, 4 = order, 5 = superfamily. Only 33
    candidates are depth 3 and 130 are depth 4, against 599 at depth 5, so a
    proportional sample would leave the shallow strata with no power. All
    depth-3 clusters are taken, depth-4 is over-sampled, and depth-5 fills the
    remainder.
  * `n_tools` — within each depth band, spread across 2/3/4/5 supporting tools
    so the depth effect is not just the support effect wearing a hat.

Selection is seeded and the stratum label is written out with each cluster, so
the batch is reproducible and the analysis can weight or report per stratum.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def path_depth(p) -> int:
    """Number of levels the canonical path resolves to ('repeat:TE:ClassI' = 3)."""
    if not isinstance(p, str) or not p:
        return 0
    return p.count(":") + 1


def select(cand: pd.DataFrame, n: int, seed: int = 42,
           depth_quota: dict[int, float] | None = None) -> pd.DataFrame:
    """Return `n` candidate clusters, stratified by path depth then n_tools."""
    df = cand[cand.candidate].copy()
    df["path_depth"] = df.majority_path.map(path_depth)
    rng = np.random.default_rng(seed)

    # take every shallow cluster; they are the scarce, informative stratum
    quota = depth_quota or {3: 1.0, 4: 0.45, 5: None}
    picked = []
    remaining = n
    for depth in (3, 4):
        pool = df[df.path_depth == depth]
        want = len(pool) if quota.get(depth) == 1.0 else int(round(
            len(pool) * float(quota.get(depth, 0.0))))
        want = min(want, len(pool), remaining)
        picked.append(_spread_by_tools(pool, want, rng))
        remaining -= want
    rest = df[~df.path_depth.isin((3, 4))]
    picked.append(_spread_by_tools(rest, min(remaining, len(rest)), rng))

    out = pd.concat(picked, ignore_index=True)
    return out.sort_values(["path_depth", "n_tools", "cluster_id"]
                           ).reset_index(drop=True)


def _spread_by_tools(pool: pd.DataFrame, want: int,
                     rng: np.random.Generator) -> pd.DataFrame:
    """Spread `want` picks evenly across n_tools, and WITHIN each n_tools level
    evenly across majority_order.

    The order pass is not cosmetic. Without it the draw follows the pool's
    composition, and the pool is dominated by short non-autonomous families:
    the first 200-cluster batch drew 2 of 19 long DNA/RC candidates against a
    stratum expectation of ~4.8, and every domain-positive DNA/RC candidate in
    the assembly (clusters 479, 755, 787, 1008) went unselected -- so the built
    batch contained no autonomous DNA transposon at all, and that read as a
    property of the genome until it was checked.
    """
    if want <= 0 or not len(pool):
        return pool.iloc[:0]
    if want >= len(pool):
        return pool
    by_tools = {k: g for k, g in pool.groupby("n_tools")}
    picks, need = [], want
    taken: dict = {}
    # round-robin over n_tools, and inside each level round-robin over order,
    # so a scarce order is never crowded out by an abundant one
    queues: dict = {}
    for k, g in by_tools.items():
        orders = {o: sub.iloc[rng.permutation(len(sub))]
                  for o, sub in g.groupby("majority_order", dropna=False)}
        queues[k] = (sorted(orders), orders, {o: 0 for o in orders})
    levels = sorted(by_tools)
    progress = True
    while need > 0 and progress:
        progress = False
        for k in levels:
            if need == 0:
                break
            names, orders, pos = queues[k]
            for o in names:
                if need == 0:
                    break
                g = orders[o]
                if pos[o] >= len(g):
                    continue
                picks.append(g.iloc[[pos[o]]])
                pos[o] += 1
                need -= 1
                progress = True
    return pd.concat(picks, ignore_index=True) if picks else pool.iloc[:0]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates",
                    default="runs/GCA_951799975.1/candidates_strict.tsv")
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="runs/GCA_951799975.1/batch_selection.tsv")
    args = ap.parse_args(argv)

    cand = pd.read_csv(args.candidates, sep="\t")
    sel = select(cand, args.n, args.seed)
    cols = [c for c in ["cluster_id", "cluster_key", "path_depth", "n_tools",
                        "n_members", "majority_order", "majority_path",
                        "pooled_full_len", "pooled_bp", "div_median"]
            if c in sel.columns]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    sel[cols].to_csv(args.out, sep="\t", index=False)

    print(f"selected {len(sel)} of {int(cand.candidate.sum())} candidates "
          f"-> {args.out}")
    print("\nby path_depth x n_tools:")
    print(pd.crosstab(sel.path_depth, sel.n_tools, margins=True).to_string())
    print("\nby majority_order:")
    print(sel.majority_order.value_counts().to_string())
    print("\npooled_full_len percentiles covered: "
          + ", ".join(f"p{p}={int(np.percentile(sel.pooled_full_len, p))}"
                      for p in (5, 25, 50, 75, 95)))


if __name__ == "__main__":
    main()

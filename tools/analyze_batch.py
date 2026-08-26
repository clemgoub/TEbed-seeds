"""Batch-level analysis of built seed packets (NEXT_INSTRUCTIONS tasks A2, A4).

    PYTHONPATH=. .venv/bin/python tools/analyze_batch.py \
        --work-dir runs/GCA_951799975.1 --out runs/GCA_951799975.1/analysis

Two questions, both of which the 24-packet batch was too small and too biased
to answer:

**A2 — is engine disagreement higher where the classification is shallower?**
F8b noticed that the two packets where MAFFT and Refiner disagreed most were
both in the one cluster whose canonical path stopped at `repeat:TE:ClassII`.
With n=24 and 6 discordant packets that is an anecdote. The test here is a
rank correlation between |len(mafft) - len(refiner)| and the depth of the
canonical path, plus a Mann-Whitney U comparing shallow (depth<=4) against
resolved (depth 5). A permutation null is used rather than a normal
approximation, so no scipy dependency and no distributional assumption.

Length disagreement is normalised by the mean of the two consensus lengths:
a 40 bp difference means something different on a 250 bp MITE and a 6 kb LTR
element, and the batch deliberately spans both.

**A4 — does the support/completeness relation hold at batch scale?**
F7 measured near-full-length rising 0.4% -> 73.4% with 1 -> 4 supporting tools,
pooled over 12 clusters that were all high-copy and mostly 4-tool. Recomputed
here over every locus in the batch, and reported per cluster as well as pooled,
because a pooled figure can be carried by a handful of huge clusters.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


# ------------------------------------------------------------------ stats
def spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Rank correlation, average ranks for ties."""
    rx = pd.Series(x).rank().to_numpy()
    ry = pd.Series(y).rank().to_numpy()
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    den = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / den) if den else float("nan")


def perm_p(x: np.ndarray, y: np.ndarray, stat, n_perm: int = 20000,
           seed: int = 0) -> tuple[float, float]:
    """Two-sided permutation p-value for `stat(x, y)`.

    Permutation rather than a table: the disagreement distribution is heavily
    zero-inflated (most packets agree exactly), so a normal approximation on
    the correlation would be optimistic.
    """
    rng = np.random.default_rng(seed)
    obs = stat(x, y)
    if not np.isfinite(obs):
        return obs, float("nan")
    null = np.empty(n_perm)
    yy = y.copy()
    for i in range(n_perm):
        rng.shuffle(yy)
        null[i] = stat(x, yy)
    p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (n_perm + 1)
    return float(obs), float(p)


def mann_whitney(a: np.ndarray, b: np.ndarray, n_perm: int = 20000,
                 seed: int = 0) -> dict:
    """U, rank-biserial effect size, and a permutation p-value."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if not len(a) or not len(b):
        return dict(n_a=len(a), n_b=len(b), U=float("nan"),
                    effect=float("nan"), p=float("nan"))
    both = np.concatenate([a, b])
    lab = np.concatenate([np.zeros(len(a)), np.ones(len(b))])

    def U_of(labels):
        r = pd.Series(both).rank().to_numpy()
        n1 = int((labels == 0).sum())
        return r[labels == 0].sum() - n1 * (n1 + 1) / 2

    obs = U_of(lab)
    rng = np.random.default_rng(seed)
    ll = lab.copy()
    null = np.empty(n_perm)
    centre = len(a) * len(b) / 2
    for i in range(n_perm):
        rng.shuffle(ll)
        null[i] = U_of(ll)
    p = (np.sum(np.abs(null - centre) >= abs(obs - centre)) + 1) / (n_perm + 1)
    return dict(n_a=len(a), n_b=len(b), U=float(obs),
                effect=float(2 * obs / (len(a) * len(b)) - 1), p=float(p))


# ------------------------------------------------------------------ loading
def load_packets(work: Path) -> pd.DataFrame:
    rows = []
    for f in sorted(work.glob("seed_packets/cluster_*/*/packet.json")):
        try:
            p = json.load(open(f))
        except Exception:
            continue
        if "error" in p or not p.get("engines"):
            continue
        eng = p["engines"]
        rec = dict(
            cluster_id=p["cluster_id"], cluster_key=p.get("cluster_key"),
            mode=p["mode"], canonical_path=p.get("canonical_path"),
            tp=p.get("tp"), tp_status=p.get("tp_status"),
            n_loci=p.get("n_loci"), n_loci_full=p.get("n_loci_full"),
            modal_cons=p.get("modal_consensus_len"),
            n_deconvolved=len(p.get("over_modal_by_member") or {}),
            support=json.dumps(p.get("support_histogram") or {}),
            full_by_support=json.dumps(p.get("full_frac_by_support") or {}),
        )
        for name in ("mafft", "refiner"):
            e = eng.get(name) or {}
            if "error" in e:
                rec[f"{name}_error"] = e["error"]
                continue
            rec[f"{name}_cons"] = e.get("rebuilt_consensus_len")
            rec[f"{name}_rows"] = e.get("alignment_rows")
            rec[f"{name}_depth"] = e.get("median_depth")
            lint = e.get("lint") or {}
            gen = lint.get("genome") or {}
            cc = lint.get("coord_check") or {}
            rec[f"{name}_lint_err"] = gen.get("n_error")
            rec[f"{name}_coord_warn"] = gen.get("n_coord_warnings")
            rec[f"{name}_verified"] = cc.get("n_verified")
            rec[f"{name}_n_rows_cc"] = cc.get("n_rows")
        rows.append(rec)
    d = pd.DataFrame(rows)
    if len(d):
        d["path_depth"] = d.canonical_path.fillna("").map(
            lambda p: p.count(":") + 1 if p else 0)
    return d


# ------------------------------------------------------------------ A2
def analyse_engine_disagreement(d: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    # A packet whose alignment has NO match column produces a zero-length
    # consensus. Counting that as "the engines disagree by 100%" is nonsense --
    # it is an engine FAILURE and belongs in a failure rate, not in a
    # disagreement distribution. Left in, it moved the median disagreement from
    # 15 bp to 82 bp and made the relative measure undefined.
    ok = d.dropna(subset=["mafft_cons", "refiner_cons"]).copy()
    ok = ok[(ok.path_depth > 0) & (ok.mafft_cons > 0) & (ok.refiner_cons > 0)]
    ok["abs_diff"] = (ok.mafft_cons - ok.refiner_cons).abs()
    ok["rel_diff"] = ok.abs_diff / ((ok.mafft_cons + ok.refiner_cons) / 2)
    ok["shallow"] = ok.path_depth <= 4

    res: dict = {"n_packets": int(len(ok)),
                 "n_exact_agree": int((ok.abs_diff == 0).sum()),
                 "median_abs_diff": float(ok.abs_diff.median()),
                 "median_rel_diff": float(ok.rel_diff.median())}
    for col in ("abs_diff", "rel_diff"):
        rho, p = perm_p(ok.path_depth.to_numpy(float), ok[col].to_numpy(float),
                        spearman)
        res[f"spearman_depth_vs_{col}"] = round(rho, 4)
        res[f"spearman_depth_vs_{col}_p"] = round(p, 5)
        mw = mann_whitney(ok.loc[ok.shallow, col].to_numpy(),
                          ok.loc[~ok.shallow, col].to_numpy())
        res[f"mw_{col}"] = {k: (round(v, 4) if isinstance(v, float) else v)
                            for k, v in mw.items()}
        res[f"median_{col}_shallow"] = float(ok.loc[ok.shallow, col].median())
        res[f"median_{col}_resolved"] = float(ok.loc[~ok.shallow, col].median())
    return ok, res


# ------------------------------------------------------------------ A4
def analyse_support(work: Path, clusters: list[int]) -> pd.DataFrame:
    """Pool loci.tsv across the batch: tool support vs near-full-length."""
    frames = []
    for cid in clusters:
        f = work / "seed_packets" / f"cluster_{cid:05d}" / "gap_aware" / "loci.tsv"
        if not f.exists():
            continue
        try:
            t = pd.read_csv(f, sep="\t", usecols=["n_tools", "any_full"])
        except Exception:
            continue
        t["cluster_id"] = cid
        frames.append(t)
    if not frames:
        return pd.DataFrame()
    all_loci = pd.concat(frames, ignore_index=True)
    pooled = (all_loci.groupby("n_tools")
              .agg(loci=("any_full", "size"), frac_full=("any_full", "mean"))
              .reset_index())
    pooled["frac_full"] = pooled.frac_full.round(4)
    # per cluster too: a pooled figure can be carried by a few huge clusters
    per = (all_loci.groupby(["cluster_id", "n_tools"]).any_full.mean()
           .reset_index().pivot(index="cluster_id", columns="n_tools",
                                values="any_full"))
    pooled.attrs["per_cluster"] = per
    return pooled


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", default="runs/GCA_951799975.1")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    work = Path(args.work_dir)
    out = Path(args.out or work / "analysis")
    out.mkdir(parents=True, exist_ok=True)

    d = load_packets(work)
    if not len(d):
        raise SystemExit(f"no packets under {work}/seed_packets")
    d.to_csv(out / "packets.tsv", sep="\t", index=False)
    print(f"packets: {len(d)}  clusters: {d.cluster_id.nunique()}")

    # health first -- an analysis of broken seeds is not worth running
    for e in ("mafft", "refiner"):
        err = d.get(f"{e}_lint_err")
        cw = d.get(f"{e}_coord_warn")
        v, n = d.get(f"{e}_verified"), d.get(f"{e}_n_rows_cc")
        if err is not None:
            print(f"  {e:8s} lint ERRORs {int(err.fillna(0).sum()):4d} | "
                  f"coord warnings {int(cw.fillna(0).sum()):4d} | "
                  f"rows verified {int(v.fillna(0).sum()):,}/"
                  f"{int(n.fillna(0).sum()):,}")

    # degenerate packets first: they are a result in their own right
    print("\n=== engines that produced NO consensus (zero match columns) ===")
    for e in ("mafft", "refiner"):
        c = d.get(f"{e}_cons")
        if c is None:
            continue
        z = int((c.fillna(0) == 0).sum())
        print(f"  {e:8s} {z:4d}/{len(d)} packets ({z/len(d):.1%})")
    deg = d[(d.mafft_cons.fillna(0) == 0) | (d.refiner_cons.fillna(0) == 0)]
    if len(deg):
        deg.to_csv(out / "degenerate_packets.tsv", sep="\t", index=False)
        print(f"  affected clusters: {deg.cluster_id.nunique()} "
              f"-> degenerate_packets.tsv")
        b = d.assign(bad=(d.mafft_cons.fillna(0) == 0)
                     | (d.refiner_cons.fillna(0) == 0))
        print("  degenerate rate by tool support (n_tools of the cluster):")
        print(b.groupby("path_depth").bad.agg(["size", "mean"]).round(3).to_string())

    ok, res = analyse_engine_disagreement(d)
    ok.to_csv(out / "engine_disagreement.tsv", sep="\t", index=False)
    json.dump(res, open(out / "a2_engine_disagreement.json", "w"), indent=1)
    print("\n=== A2: engine disagreement vs classification depth ===")
    print(f"  packets compared        {res['n_packets']}")
    print(f"  exact agreement         {res['n_exact_agree']} "
          f"({res['n_exact_agree'] / max(res['n_packets'],1):.1%})")
    print(f"  median |diff|           {res['median_abs_diff']:.1f} bp "
          f"({res['median_rel_diff']:.3%} of consensus length)")
    for col in ("abs_diff", "rel_diff"):
        mw = res[f"mw_{col}"]
        print(f"  --- {col} ---")
        print(f"    Spearman(depth, {col}) = {res[f'spearman_depth_vs_{col}']:+.3f}"
              f"   permutation p = {res[f'spearman_depth_vs_{col}_p']:.4f}")
        print(f"    shallow (depth<=4) median {res[f'median_{col}_shallow']:.4g}"
              f"  vs resolved (5) {res[f'median_{col}_resolved']:.4g}")
        print(f"    Mann-Whitney U={mw['U']:.0f} n={mw['n_a']}/{mw['n_b']} "
              f"effect={mw['effect']:+.3f} p={mw['p']:.4f}")

    print("\n  by path depth:")
    print(ok.groupby("path_depth").agg(
        packets=("abs_diff", "size"), median_abs=("abs_diff", "median"),
        median_rel=("rel_diff", "median"),
        frac_exact=("abs_diff", lambda s: (s == 0).mean())).round(4).to_string())

    pooled = analyse_support(work, sorted(d.cluster_id.unique()))
    if len(pooled):
        pooled.to_csv(out / "a4_support_vs_completeness.tsv", sep="\t", index=False)
        print("\n=== A4: tool support vs near-full-length, pooled over the batch ===")
        print(pooled.to_string(index=False))
        per = pooled.attrs["per_cluster"]
        per.to_csv(out / "a4_support_per_cluster.tsv", sep="\t")
        mono = per.apply(lambda r: r.dropna().is_monotonic_increasing, axis=1)
        print(f"\n  clusters whose frac_full rises monotonically with n_tools: "
              f"{int(mono.sum())}/{len(mono)}")
    print(f"\nwritten to {out}")


if __name__ == "__main__":
    main()

"""Stage 1 -- copies: gap-aware merge, extraction windows, sampling.

Coordinate-only: everything here runs from the BED16s; sequence extraction
happens later against the assembly FASTA (config: assembly_fasta).

Measured basis (PLAN_A rev 3 section 0.3): 62.9% of rm2 inter-fragment gaps
are mostly occupied by ANOTHER element (median gap 2.2 kb) -- merge-always
would absorb nested insertions into the copy. Rule: merge across a gap only
when the gap is unclaimed, or claimed by a family in the SAME cluster.
Both modes are implemented so the first batch can compare them.
"""
from __future__ import annotations

import collections

import numpy as np
import pandas as pd

MAX_MERGE_GAP = 20_000     # beyond this, never merge (element-scale sanity cap)
GAP_FOREIGN_FRAC = 0.5     # gap "claimed" when >50% covered by out-of-cluster hits
MAX_COPY_X_CONSENSUS = 1.5  # a copy longer than this x consensus is not one copy


def consensus_length(hits: pd.DataFrame) -> float:
    """MODAL implied consensus length (repeat_end + repeat_left).

    NOT max: families flagged len_inconsistent upstream carry a handful of
    outlier rows whose implied length is an order of magnitude too large
    (measured: rnd-1_family-286 -> modal 269 bp from 31,387 hits vs max
    4,202 bp from 11 hits). Using max made near-full-length counts read ZERO
    where the true count is 14,642.
    """
    implied = (hits.repeat_end + hits.repeat_left).dropna()
    if not len(implied):
        return float("nan")
    return float(implied.mode().iat[0])


def build_foreign_probe(repo, tools, chrom_sizes, cluster_fams, chroms=None):
    """Pre-paint per-chrom boolean coverage of OUT-OF-CLUSTER hits, once.

    The naive alternative (re-scan every tool's BED per gap) is what made the
    first stage-1 run intractable: ~6M rows x tens of thousands of gaps.
    Painting once per chromosome makes each gap check an array-slice mean.
    Returns (cov_dict, probe_fn) where probe_fn(chrom, a, b) -> foreign frac.
    """
    from .stage0 import BED16_COLS
    want = set(chroms) if chroms else None
    cov: dict[str, np.ndarray] = {}
    for t in tools:
        df = pd.read_csv(repo / "inputs" / f"{t}.bed", sep="\t", names=BED16_COLS,
                         skiprows=1, usecols=["chrom", "chromStart", "chromEnd", "name"])
        if want:
            df = df[df.chrom.isin(want)]
        own = {f for (tt, f) in cluster_fams if tt == t}
        df = df[~df.name.isin(own)]           # same-cluster families are not foreign
        for ch, g in df.groupby("chrom"):
            L = chrom_sizes.get(ch)
            if L is None:
                continue
            a = cov.get(ch)
            if a is None:
                a = np.zeros(L, dtype=bool)
                cov[ch] = a
            for s, e in zip(g.chromStart.to_numpy(), g.chromEnd.to_numpy()):
                a[s:e] = True

    def probe(ch, a, b):
        arr = cov.get(ch)
        if arr is None or b <= a:
            return 0.0
        return float(arr[a:b].mean())

    return cov, probe


def build_clustermate_conslen(repo, cluster_fams, chrom_sizes,
                              conslen_by_member: dict[str, float],
                              weight_by_member: dict[str, float] | None = None):
    """Paint cluster-mates' consensus length per base, for inheritance.

    A member without consensus coordinates (EDTA: 0% of hits carry them)
    cannot judge near-full-length on its own. PLAN_A 4.1 resolves this
    per LOCUS, not per cluster -- and that distinction matters: in the
    measured cluster the mates' consensus lengths are bimodal (264/269 vs
    764/765 bp), so any cluster-wide average would be wrong for both groups.
    Returns lookup(chrom, start, end) -> modal cluster-mate consensus length
    over that span, or nan.

    ORDER MATTERS AND MUST NOT BE INCIDENTAL. Members are painted into a shared
    array, so at a base two mates both cover, whoever paints LAST wins. Callers
    used to pass a set, whose iteration order depends on PYTHONHASHSEED -- so
    the inherited length at contested bases changed between runs on identical
    input. Measured: one EDTA member's inherited consensus length moved from
    381 bp to 1,755 bp between two runs with no other change. That is the
    cluster_id reproducibility bug again, one level down, and this time it
    moves a scientific value rather than a label.

    Members are now painted in order of INCREASING evidence (hit count by
    default), so the best-supported mate paints last and wins a contested base.
    """
    order = sorted(cluster_fams,
                   key=lambda tf: (weight_by_member or {}).get(
                       f"{tf[0]}:{tf[1]}", 0.0))
    from .stage0 import BED16_COLS
    arrs: dict[str, np.ndarray] = {}
    for (tool, fam) in order:
        L = conslen_by_member.get(f"{tool}:{fam}")
        if L is None or np.isnan(L):
            continue
        df = pd.read_csv(repo / "inputs" / f"{tool}.bed", sep="\t", names=BED16_COLS,
                         skiprows=1, usecols=["chrom", "chromStart", "chromEnd", "name"])
        df = df[df.name == fam]
        for ch, g in df.groupby("chrom"):
            n = chrom_sizes.get(ch)
            if n is None:
                continue
            a = arrs.get(ch)
            if a is None:
                a = np.zeros(n, dtype=np.int32)
                arrs[ch] = a
            for s, e in zip(g.chromStart.to_numpy(), g.chromEnd.to_numpy()):
                a[s:e] = int(L)

    def lookup(ch, s, e):
        a = arrs.get(ch)
        if a is None:
            return float("nan")
        seg = a[s:e]
        seg = seg[seg > 0]
        if not seg.size:
            return float("nan")
        vals, cts = np.unique(seg, return_counts=True)
        return float(vals[cts.argmax()])

    return lookup


def member_hits(repo, tool: str, family: str) -> pd.DataFrame:
    """All hits of one member family, with consensus coords where present."""
    from .stage0 import BED16_COLS
    df = pd.read_csv(
        repo / "inputs" / f"{tool}.bed", sep="\t", names=BED16_COLS, skiprows=1,
        usecols=["chrom", "chromStart", "chromEnd", "name", "strand",
                 "perc_div", "repeat_start", "repeat_end", "repeat_left", "hit_id"],
        na_values=["NA"])
    return df[df.name == family].reset_index(drop=True)


def merge_copies(hits: pd.DataFrame, mode: str,
                 foreign_cov=None, cons_len: float = float("nan")) -> pd.DataFrame:
    """Group hits into copies.

    mode='hit_id'      -- trust the tool's own fragment linkage (REPET).
    mode='merge_always'-- adjacent same-family hits merge across any gap <= cap.
    mode='gap_aware'   -- merge only if the gap is NOT mostly claimed by
                          out-of-cluster sequence (foreign_cov callable:
                          (chrom, a, b) -> fraction of [a,b) covered).
    Emits one row per copy: chrom, start, end, strand, n_fragments,
    frag_starts/frag_ends (';'-joined -- extraction uses fragments, never the
    gap), plus consensus-coord span where available.
    """
    rows = []
    if mode == "hit_id" and hits.hit_id.notna().any():
        groups = hits.groupby("hit_id", sort=False)
    else:
        # positional grouping per chrom+strand
        h = hits.sort_values(["chrom", "chromStart"]).reset_index(drop=True)
        gid = np.zeros(len(h), dtype=int)
        g = 0
        for i in range(1, len(h)):
            same = (h.chrom[i] == h.chrom[i - 1]) and (h.strand[i] == h.strand[i - 1])
            gap_a, gap_b = h.chromEnd[i - 1], h.chromStart[i]
            gap = gap_b - gap_a
            merge = same and 0 <= gap <= MAX_MERGE_GAP
            if merge and mode == "gap_aware" and gap > 0 and foreign_cov is not None:
                merge = foreign_cov(h.chrom[i], gap_a, gap_b) <= GAP_FOREIGN_FRAC
            if not merge:
                g += 1
            gid[i] = g
        h["_g"] = gid
        # Runaway-chain guard: transitive merging can walk a whole tandem array
        # of an element into one "copy" (measured: 328x the consensus length,
        # 32% of merge_always copies over 1.5x). A copy cannot exceed the
        # consensus by more than MAX_COPY_X_CONSENSUS -- split the group at its
        # widest internal gap until every piece is within bound.
        if not np.isnan(cons_len):
            cap = MAX_COPY_X_CONSENSUS * cons_len
            next_g = int(h["_g"].max()) + 1
            changed = True
            while changed:
                changed = False
                for gv, grp in list(h.groupby("_g", sort=False)):
                    if len(grp) < 2:
                        continue
                    if grp.chromEnd.max() - grp.chromStart.min() <= cap:
                        continue
                    gi = grp.sort_values("chromStart")
                    gaps = gi.chromStart.to_numpy()[1:] - gi.chromEnd.to_numpy()[:-1]
                    cut = int(np.argmax(gaps)) + 1
                    h.loc[gi.index[cut:], "_g"] = next_g
                    next_g += 1
                    changed = True
        groups = h.groupby("_g", sort=False)

    for _, grp in groups:
        grp = grp.sort_values("chromStart")
        rs = grp.repeat_start.min() if grp.repeat_start.notna().any() else np.nan
        re_ = grp.repeat_end.max() if grp.repeat_end.notna().any() else np.nan
        rows.append(dict(
            chrom=grp.chrom.iat[0], start=int(grp.chromStart.min()),
            end=int(grp.chromEnd.max()), strand=grp.strand.iat[0],
            n_fragments=len(grp),
            frag_starts=";".join(map(str, grp.chromStart)),
            frag_ends=";".join(map(str, grp.chromEnd)),
            div=float(np.nanmean(grp.perc_div)) if grp.perc_div.notna().any() else np.nan,
            cons_start=rs, cons_end=re_,
        ))
    return pd.DataFrame(rows)


def near_full_length(copies: pd.DataFrame, cons_len: float,
                     span_frac=0.8, span_only: bool = False) -> pd.Series:
    """Both consensus ends reached (within 20 bp) and >=span_frac of length.

    `span_only` drops the end-reaching test and asks only that the copy spans
    span_frac of a unit. It is for members whose consensus was DECONVOLVED as
    an over-assembly: their repeat_start/repeat_end refer to the multi-unit
    coordinate system of the collapsed entry, so "reached both ends" is not a
    meaningful question about a single unit, and asking it against the unit
    length would pass almost everything.

    Measured on the four confirmed over-assembled REPET families (10,229 hits),
    against a cross-tool ground truth of 49.5% near-full-length taken from the
    independent rm2 consensus of the same element:
        own 764-771 bp consensus, ends test    0.74%   (the old behaviour, ~67x low)
        unit length, ends test kept           25.2%    (wrong: mixed coordinate systems)
        span-only against the unit length     58.2%    (closest to truth)
    """
    if np.isnan(cons_len) or copies.cons_start.isna().all():
        return pd.Series(False, index=copies.index)
    span = (copies.cons_end - copies.cons_start) / cons_len
    if span_only:
        return span >= span_frac
    return (copies.cons_start <= 20) & (copies.cons_end >= cons_len - 20) & (span >= span_frac)


def sample_copies(copies: pd.DataFrame, is_full: pd.Series, cap: int,
                  floor: int, full_frac: float, seed=42) -> pd.DataFrame:
    """Stratified sampling: ~full_frac near-full-length, rest spanning copies.

    Realized (not target) composition is what gets recorded -- repeatome
    intactness varies by taxon, so the achievable ratio does too.
    """
    rng = np.random.default_rng(seed)
    if len(copies) <= max(floor, 1):
        out = copies.copy()
        out["sampled_as"] = np.where(is_full, "full", "partial")
        return out
    n = min(cap, len(copies))
    want_full = int(round(n * full_frac))
    fl = copies[is_full]
    pt = copies[~is_full]
    take_full = fl.sample(min(want_full, len(fl)), random_state=rng.integers(1 << 31))
    rest = n - len(take_full)
    # prefer partials that extend coverage toward the consensus edges
    if len(pt) and rest > 0:
        take_part = pt.sample(min(rest, len(pt)), random_state=rng.integers(1 << 31))
    else:
        take_part = pt.iloc[:0]
    short = n - len(take_full) - len(take_part)
    if short > 0 and len(fl) > len(take_full):
        extra = fl.drop(take_full.index).sample(
            min(short, len(fl) - len(take_full)), random_state=rng.integers(1 << 31))
        take_full = pd.concat([take_full, extra])
    out = pd.concat([take_full.assign(sampled_as="full"),
                     take_part.assign(sampled_as="partial")])
    return out.sort_values(["chrom", "start"]).reset_index(drop=True)


def smitten(assembly: str, row) -> str:
    """Dfam Smitten id: 1-based fully-closed, from BED half-open."""
    return f"{assembly}:{row.chrom}:{row.start + 1}-{row.end}_{row.strand}"

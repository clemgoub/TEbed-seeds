"""Cluster-level candidate gates (PLAN_A rev 3 section 3.4).

Each gate is reported individually; the shortlist requires all. Per-locus
class conflict is never a filter -- G3 asks the cluster-level question (do the
tools broadly agree what KIND of element this is), members that stop above
order level ABSTAIN, and evidence weighting means a fragmentary member cannot
outvote the bulk.
"""
from __future__ import annotations

import collections

import numpy as np
import pandas as pd


def order_coherence(members: list[str], fam_order: dict[str, str],
                    weights: dict[str, float]) -> tuple[float, str | None]:
    """Member-weighted agreement at order level. Abstainers excluded from
    the denominator. Returns (coherence, majority_order)."""
    votes: dict[str, float] = collections.defaultdict(float)
    for m in members:
        o = fam_order.get(m)
        if o in (None, "repeat"):     # stops above order: abstains
            continue
        votes[o] += weights.get(m, 1.0)
    if not votes:
        return 0.0, None
    tot = sum(votes.values())
    best = max(votes, key=votes.get)
    return votes[best] / tot, best


def majority_path(members: list[str], fam_path: dict[str, str],
                  weights: dict[str, float]) -> tuple[str | None, float]:
    """Member-weighted majority canonical classification path.

    Carried through to stage 2 because `#=GF TP` needs a classification and
    `majority_order` alone ("LINE") is too coarse to map onto the Dfam
    vocabulary. Uninformative paths abstain, exactly as order coherence does,
    so a member that stops at "repeat" cannot outvote one that resolved a
    superfamily.
    """
    votes: dict[str, float] = collections.defaultdict(float)
    for m in members:
        p = fam_path.get(m)
        if not p or p in ("repeat", "repeat:TE"):
            continue
        votes[p] += weights.get(m, 1.0)
    if not votes:
        return None, 0.0
    best = max(votes, key=votes.get)
    return best, votes[best] / sum(votes.values())


def evaluate_clusters(clusters: list[list[str]], evidence: pd.DataFrame,
                      fam_order: dict[str, str], cfg: dict,
                      fam_path: dict[str, str] | None = None) -> pd.DataFrame:
    """evidence: one row per 'tool:family' member with columns
    n_full_len, genomic_bp, cov_ge3_frac, div_median, frac_tandemtool.
    Missing members (e.g. EDTA families absent from the upstream evidence
    table) contribute copies but not full-length counts -- they are NaN-safe.
    """
    ev = evidence.set_index("member")
    g = cfg["gates"]
    rows = []
    for ci, members in enumerate(clusters):
        sub = ev.reindex(members)
        w = dict(zip(members, sub.genomic_bp.fillna(0.0) + 1.0))
        coh, majority = order_coherence(members, fam_order, w)
        mpath, mpath_frac = majority_path(members, fam_path or {}, w)
        pooled_full = float(np.nansum(sub.n_full_len))
        pooled_bp = float(np.nansum(sub.genomic_bp))
        tandem_frac = (float(np.nansum(sub.frac_tandemtool * sub.genomic_bp)) / pooled_bp
                       if pooled_bp else float("nan"))
        div_med = float(np.nanmedian(sub.div_median)) if sub.div_median.notna().any() else float("nan")
        depth_ok = bool((sub.cov_ge3_frac >= 0.999).any())  # any member spans at depth>=3; pooled depth refined in stage 2
        n_tools = len({m.split(":")[0] for m in members})
        gates = dict(
            G1_full_length=pooled_full >= g["min_pooled_full_length"],
            G2_depth=depth_ok,
            G3_order_coherence=coh >= cfg["order_coherence_min"],
            G4_tandem=(not np.isnan(tandem_frac)) and tandem_frac <= g["max_tandem_overlap"],
            G5_divergence=(not np.isnan(div_med)) and div_med <= g["max_median_divergence"],
            G6_tools=n_tools >= g["min_tools"],
        )
        rows.append(dict(
            cluster_id=ci, n_members=len(members), n_tools=n_tools,
            members=";".join(members), majority_order=majority,
            majority_path=mpath, majority_path_frac=round(mpath_frac, 3),
            order_coherence=round(coh, 3), pooled_full_len=int(pooled_full),
            pooled_bp=int(pooled_bp), tandem_frac=round(tandem_frac, 3) if not np.isnan(tandem_frac) else np.nan,
            div_median=round(div_med, 2) if not np.isnan(div_med) else np.nan,
            **gates, n_gates=sum(gates.values()), candidate=all(gates.values()),
        ))
    return pd.DataFrame(rows)

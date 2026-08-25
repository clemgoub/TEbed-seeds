#!/usr/bin/env python3
"""Sequence-free detection of OVER-ASSEMBLED library consensi.

An over-assembled consensus is one library entry that is really several tandem
copies of a shorter element collapsed together (RESUME.md section 5,
PIPELINE_FINDINGS.md F4).  Treating its advertised length as the element
length makes every genuine copy look ~(1 - 1/k) truncated and sinks the family
at any near-full-length gate.

The detector needs NO library FASTA.  It uses only BED16 columns:

    implied consensus length  L = repeat_end + repeat_left   (MODAL, not max
                                  -- stage1.consensus_length / F6)
    consensus match interval  [repeat_start, repeat_end]
    genomic hit length        chromEnd - chromStart

=============================================================================
THE RULE
=============================================================================
A family is flagged for a small integer k in K_RANGE (default 2..5) when ALL
THREE of the following hold.  (a) and (b) are the two conditions briefed in
RESUME.md section 5 task 4; (c) was added after (a)+(b) alone were measured to
have no specificity -- see CALIBRATION below.

(a) LENGTH COMMENSURABILITY.  The median genomic hit length m satisfies

        RHO_LO <= rho_k <= RHO_HI ,      rho_k = m * k / L

    i.e. the typical copy is one k-th of the advertised consensus.
    Defaults RHO_LO = 0.70, RHO_HI = 1.30.  Justification:
      * upper bound 1.30 -- if copies really are single units they cannot be
        much LONGER than a unit; 30% headroom absorbs indels, alignment slop
        and the fact that the modal L can under-shoot the true entry length.
      * lower bound 0.70 -- the measured positives (REPET, cluster 62) sit at
        rho_3 = 0.884 and 0.871.  A 0.70 floor still admits families ~20% more
        degraded than the measured instance while excluding families whose
        copies are sub-unit fragments; for those, m no longer identifies k
        (rm2 rnd-5_family-3676: L = 2351, m = 88, rho_5 = 0.19).
      * the same window puts the rm2 and pantera consensi of the SAME element
        (269 / 264 bp) at rho_3 = 2.50 / 2.69 -- outside by ~2x.
    Alone, (a) is weak: 63.6% of testable families satisfy it for some k in
    2..5, because every degraded family has a small m/L.

(b) BLOCK STRUCTURE, not uniform spread.
    Statistic = PHASE-0 ENRICHMENT OF INTERNAL MATCH ENDPOINTS.

      1. Endpoint pool E: every repeat_start and every repeat_end lying
         strictly INSIDE the consensus, in (tau, L - tau) with
         tau = max(TAU_ABS, TAU_FRAC * L)  (defaults 20 bp, 5%).
         Terminal endpoints are dropped: every family piles endpoints at 0 and
         at L, those carry no information about internal structure, and they
         would sit at phase 0 for every integer k by construction.
      2. For period P = L/k, phase phi = (x mod P)/P, binned into NB = 10 bins
         with bin 0 CENTRED on phase 0.  All k-1 internal block boundaries
         j*L/k fold onto phase 0.
      3. phase_enrich_k = h[0] / (n/NB)
         phase_z_k      = (h[0] - n/NB) / sqrt(n * (1/NB) * (1 - 1/NB))

    NULL.  Under "internal endpoints carry no k-periodic structure" the phase
    is uniform on [0,1), so h[0] ~ Binomial(n, 1/NB), E[phase_enrich] = 1 and
    phase_z is standard normal.  That parametric null is ANTI-conservative in
    practice (real endpoint distributions are lumpy for reasons unrelated to
    k), so the threshold is placed on the EFFECT SIZE:

        (b) holds iff  phase_enrich_k >= ENRICH_MIN (1.5)
                  and  phase_z_k      >= Z_MIN (5.0)

    The z floor only guards small n; ENRICH_MIN does the discrimination.
    Bin 0 is ANCHORED on phase 0 rather than taken as the maximum bin on
    purpose: an offset-free "strongest phase bin" version fires on any family
    with a single internal endpoint hotspot anywhere (measured: rm2
    rnd-1_family-286, a genuine 269 bp TIR element, reaches 2.57 that way at
    k=2).  Anchoring demands the hotspots sit at the consensus-relative
    positions j*L/k, which is what "the entry is k units" means.
    Alone, (b) holds for 16.1% of testable families for some k in 2..5.

(c) THE ADVERTISED LENGTH IS NEVER REALISED.

        full_span_frac = fraction of hits with
                         repeat_start <= 20  AND  repeat_end >= L - 20
                         AND (repeat_end - repeat_start)/L >= 0.80
                       < FULL_SPAN_MAX  (default 0.02)

    This is exactly stage1.near_full_length evaluated on the raw hits.  If a
    family HAS near-full-length copies, its consensus is realised somewhere in
    the genome and it is not an over-assembly, however block-structured its
    partial hits look.  This is the condition that ties the flag directly to
    the pipeline consequence: a flagged family is precisely one whose own
    cons_len can never produce a near-full-length copy.
    Measured: rm2 rnd-1_family-286 full_span_frac = 0.494, pantera
    Unknown_572-fGobNig = 0.573 (both excluded); the two REPET positives
    0.0055 and 0.0067 (both retained).

If several k satisfy all three, the reported k is the one with the largest
phase_enrich.  unit_len = L / k is the deconvolved element length.

=============================================================================
CALIBRATION -- and an honest statement of what this can and cannot decide
=============================================================================
Specificity was measured with a DECOY NULL: the identical criterion evaluated
at HALF-INTEGER k' in {2.5, 3.5, 4.5}.  No tandem over-assembly can produce a
non-integer number of units, so every decoy hit is a false positive of the
"the consensus is k tandem units" claim, with the family's own endpoint
distribution, copy-length distribution and sample size held fixed.

Measured on 3,650 testable families of GCA_951799975.1 (rm2, pantera,
fastltr, repet):

    criterion                       real k in 2..5      decoy k' in 2.5/3.5/4.5
    (a) only                        2322  (63.6%)       1975  (54.1%)
    (b) only                         589  (16.1%)        752  (20.6%)
    (a) and (b)                      151  ( 4.1%)        163  ( 4.5%)
    (a) and (b) and (c)               69  ( 1.9%)         55  ( 1.5%)

So (a)+(b) as briefed has NO specificity for integer k -- the decoy fires
slightly more often than the real thing.  Adding (c) produces the first real,
if modest, excess (69 vs 55; implied FDR ceiling ~80% for the integer-k
claim).

WHAT THE FLAG SUPPORTS: "this consensus is block-structured, no copy realises
its advertised length, and the length actually realised is about L/k."  That
is enough for the pipeline consequence (do not use cons_len for the
near-full-length judgement).
WHAT IT DOES NOT SUPPORT ON ITS OWN: that k is exactly an integer, i.e. that
the entry is k tandem units rather than, say, a chimera of two unrelated
elements joined at a position that happens to be commensurate with the copy
length.  Distinguishing those needs the library SEQUENCE (self-alignment:
tools/overassembly_seqtest.sh).  Downstream code should therefore prefer the
MEASURED unit length (cluster-mate inheritance) over dividing by the detected
k -- see the Part C recommendation.

boundary_support reports how many of the k-1 internal boundaries carry >= 5%
of internal endpoints.  Both known positives are single-boundary (each of the
two redundant REPET entries populates a different one of the two boundaries of
the same 3-unit array), so multi-boundary support is a stricter claim, not a
prerequisite.

=============================================================================
GUARDS
=============================================================================
  * MIN_HITS (50) hits carrying consensus coordinates.
  * MIN_CONSLEN (100 bp): below this L/k is under 20-50 bp and coordinate
    rounding dominates the phase.
  * conslen_frac >= CONSLEN_FRAC_MIN (0.5): fraction of hits whose implied
    length is within 1% of the modal L.  A family whose implied length is not
    self-consistent has no trustworthy L to divide.
  * A tool with no consensus coordinates at all (EDTA: 0 / 1,394,534 hits) is
    reported as UNTESTABLE, never as zero flagged.

=============================================================================
USAGE
=============================================================================
    python tools/detect_overassembly.py \
        --bed-dir ~/Documents/VGP_TEbed/inputs \
        --tools rm2,edta,pantera,fastltr,repet \
        --out runs/GCA_951799975.1/overassembly.tsv \
        --cache-dir /tmp/oa_cache            # optional parquet column cache
        --decoy                              # also run the half-integer null

    # importable
    from detect_overassembly import scan_tool, family_stats, PARAMS
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- parameters
BED16_COLS = [
    "chrom", "chromStart", "chromEnd", "name", "score", "strand", "SW_score",
    "perc_div", "perc_del", "perc_ins", "query_left", "repeat_class_family",
    "repeat_start", "repeat_end", "repeat_left", "hit_id",
]
USECOLS = ["chromStart", "chromEnd", "name", "repeat_class_family",
           "repeat_start", "repeat_end", "repeat_left"]
KEEPCOLS = ["name", "repeat_class_family", "repeat_start", "repeat_end",
            "repeat_left", "hit_len"]

DECOY_K = (2.5, 3.5, 4.5)   # non-integer -> every hit is a false positive


@dataclass(frozen=True)
class Params:
    k_range: tuple = (2, 3, 4, 5)
    n_bins: int = 10            # phase bins per period; bin 0 centred on 0
    tau_abs: float = 20.0       # terminal exclusion, bp
    tau_frac: float = 0.05      # terminal exclusion, fraction of L
    rho_lo: float = 0.70        # (a) window on median_len * k / L
    rho_hi: float = 1.30
    enrich_min: float = 1.5     # (b) effect size
    z_min: float = 5.0          # (b) significance floor
    full_span_max: float = 0.02  # (c) max fraction of near-full-length hits
    full_end_slop: float = 20.0  # (c) stage1.near_full_length end tolerance
    full_span_frac: float = 0.80  # (c) stage1.near_full_length span fraction
    occ_min: float = 0.05       # boundary counted as occupied at >= this share
    min_hits: int = 50
    min_conslen: float = 100.0
    conslen_frac_min: float = 0.50


PARAMS = Params()

REJECT_REASONS = ("too_few_hits_with_cons_coords", "consensus_too_short",
                  "consensus_length_inconsistent", "too_few_internal_endpoints")


# --------------------------------------------------------------------------- core stats
def modal_consensus_length(repeat_end: np.ndarray,
                           repeat_left: np.ndarray) -> tuple[float, float]:
    """MODAL implied consensus length, and the share of hits supporting it.

    Modal, not max: a handful of outlier rows imply lengths an order of
    magnitude too large -- rm2 rnd-1_family-286: modal 269 bp from 31,387 hits
    vs max 4,202 bp from 11 hits (PIPELINE_FINDINGS F6).
    """
    implied = repeat_end + repeat_left
    implied = implied[np.isfinite(implied)]
    if implied.size == 0:
        return float("nan"), float("nan")
    vals, cts = np.unique(np.round(implied).astype(np.int64), return_counts=True)
    L = float(vals[cts.argmax()])
    if L <= 0:
        return float("nan"), float("nan")
    frac = float((np.abs(implied - L) <= max(1.0, 0.01 * L)).mean())
    return L, frac


def internal_endpoints(repeat_start: np.ndarray, repeat_end: np.ndarray,
                       L: float, p: Params = PARAMS) -> np.ndarray:
    """Match endpoints strictly inside the consensus (termini carry no info)."""
    tau = max(p.tau_abs, p.tau_frac * L)
    pool = np.concatenate([repeat_start, repeat_end])
    return pool[(pool > tau) & (pool < L - tau)]


def phase_enrichment(endpoints: np.ndarray, L: float, k: float,
                     n_bins: int = PARAMS.n_bins
                     ) -> tuple[float, float, int, np.ndarray]:
    """Enrichment of internal endpoints in the phase bin centred on j*L/k.

    Null: phi ~ U[0,1) so the boundary bin holds Binomial(n, 1/n_bins).
    Accepts non-integer k (the decoy null).
    """
    n = int(endpoints.size)
    if n == 0 or not np.isfinite(L) or L <= 0:
        return float("nan"), float("nan"), 0, np.zeros(n_bins, dtype=int)
    P = L / k
    phi = ((endpoints % P) / P + 0.5 / n_bins) % 1.0     # bin 0 centred on 0
    h = np.bincount((phi * n_bins).astype(int), minlength=n_bins)[:n_bins]
    exp = n / n_bins
    return (float(h[0] / exp),
            float((h[0] - exp) / np.sqrt(exp * (1.0 - 1.0 / n_bins))), n, h)


def boundary_occupancy(endpoints: np.ndarray, L: float, k: int,
                       half_width_frac: float = 0.05) -> list[float]:
    """Share of internal endpoints within +-half_width_frac*P of each j*L/k."""
    P = L / k
    w = half_width_frac * P
    n = max(int(endpoints.size), 1)
    return [round(float((np.abs(endpoints - j * P) <= w).sum()) / n, 3)
            for j in range(1, int(k))]


def full_span_fraction(repeat_start: np.ndarray, repeat_end: np.ndarray,
                       L: float, p: Params = PARAMS) -> float:
    """Share of hits that are near-full-length (stage1.near_full_length)."""
    if not np.isfinite(L) or L <= 0 or repeat_start.size == 0:
        return float("nan")
    return float((((repeat_start <= p.full_end_slop) &
                   (repeat_end >= L - p.full_end_slop) &
                   (((repeat_end - repeat_start) / L) >= p.full_span_frac))).mean())


def family_stats(hits: pd.DataFrame, p: Params = PARAMS,
                 k_range=None) -> dict:
    """Full per-family record.  `hits` needs repeat_start/end/left + hit_len.

    Pass k_range=DECOY_K to evaluate the half-integer decoy null instead.
    """
    ks = p.k_range if k_range is None else k_range
    rs = hits["repeat_start"].to_numpy(dtype=float)
    re_ = hits["repeat_end"].to_numpy(dtype=float)
    rl = hits["repeat_left"].to_numpy(dtype=float)
    hl = hits["hit_len"].to_numpy(dtype=float)
    ok = np.isfinite(rs) & np.isfinite(re_) & np.isfinite(rl)
    n_cc = int(ok.sum())

    rec = dict(n_hits=int(len(hits)), n_hits_conscoord=n_cc,
               cons_len=float("nan"), conslen_frac=float("nan"),
               median_len=float(np.median(hl)) if hl.size else float("nan"),
               p90_len=float(np.quantile(hl, 0.90)) if hl.size else float("nan"),
               full_span_frac=float("nan"), n_internal_endpoints=0,
               k=0, unit_len=float("nan"), rho=float("nan"),
               phase_enrich=float("nan"), phase_z=float("nan"),
               boundary_occ="", boundary_support=0,
               best_enrich_any_k=float("nan"), best_enrich_k=0,
               flag=False, reason="")

    if n_cc < p.min_hits:
        rec["reason"] = "too_few_hits_with_cons_coords"
        return rec

    rs, re_, rl, hl = rs[ok], re_[ok], rl[ok], hl[ok]
    L, cfrac = modal_consensus_length(re_, rl)
    rec["cons_len"], rec["conslen_frac"] = L, cfrac
    rec["median_len"] = float(np.median(hl))
    rec["p90_len"] = float(np.quantile(hl, 0.90))

    if not np.isfinite(L) or L < p.min_conslen:
        rec["reason"] = "consensus_too_short"
        return rec
    if cfrac < p.conslen_frac_min:
        rec["reason"] = "consensus_length_inconsistent"
        return rec

    rec["full_span_frac"] = round(full_span_fraction(rs, re_, L, p), 4)
    E = internal_endpoints(rs, re_, L, p)
    rec["n_internal_endpoints"] = int(E.size)
    if E.size < p.min_hits:
        rec["reason"] = "too_few_internal_endpoints"
        return rec

    m = float(np.median(hl))
    cond_c = rec["full_span_frac"] < p.full_span_max
    best_any, passing = (-np.inf, 0), []
    for k in ks:
        enrich, z, _, _ = phase_enrichment(E, L, k, p.n_bins)
        if enrich > best_any[0]:
            best_any = (enrich, k)
        rho = m * k / L
        cond_a = p.rho_lo <= rho <= p.rho_hi
        cond_b = (enrich >= p.enrich_min) and (z >= p.z_min)
        if cond_a and cond_b and cond_c:
            passing.append((enrich, k, rho, z))
    rec["best_enrich_any_k"] = round(float(best_any[0]), 3)
    rec["best_enrich_k"] = float(best_any[1])

    if not passing:
        rec["reason"] = ("consensus_length_is_realised" if not cond_c
                         else "no_k_passes")
        return rec

    enrich, k, rho, z = max(passing)
    occ = (boundary_occupancy(E, L, int(k)) if float(k).is_integer() else [])
    rec.update(k=k, unit_len=round(L / k, 1), rho=round(float(rho), 3),
               phase_enrich=round(float(enrich), 3), phase_z=round(float(z), 2),
               boundary_occ=";".join(str(v) for v in occ),
               boundary_support=int(sum(1 for o in occ if o >= p.occ_min)),
               flag=True, reason="over_assembly")
    return rec


# --------------------------------------------------------------------------- driver
def load_tool(bed: Path, cache: Path | None = None) -> pd.DataFrame:
    """Read only the BED16 columns needed (the BEDs are 120-240 MB each)."""
    if cache is not None and cache.exists():
        df = pd.read_parquet(cache)
        if "hit_len" not in df.columns:
            df["hit_len"] = df["chromEnd"] - df["chromStart"]
        return df[[c for c in KEEPCOLS if c in df.columns]]
    parts = []
    for ch in pd.read_csv(bed, sep="\t", names=BED16_COLS, skiprows=1,
                          usecols=USECOLS, na_values=["NA"], chunksize=2_000_000,
                          dtype={"name": "string", "repeat_class_family": "string"}):
        for c in ("repeat_start", "repeat_end", "repeat_left"):
            ch[c] = pd.to_numeric(ch[c], errors="coerce")
        ch["hit_len"] = ch["chromEnd"] - ch["chromStart"]
        parts.append(ch[KEEPCOLS])
    df = pd.concat(parts, ignore_index=True)
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(cache, index=False)
    return df


def scan_tool(tool: str, bed: Path, p: Params = PARAMS,
              cache: Path | None = None, eligible: set[str] | None = None,
              k_range=None) -> pd.DataFrame:
    """Per-family records for one tool BED.  frac_cons_coords in .attrs."""
    df = load_tool(bed, cache)
    df["name"] = df["name"].astype(str)
    if eligible is not None:
        df = df[df["name"].isin(eligible)]
    has_cc = df["repeat_end"].notna() & df["repeat_left"].notna()
    frac_cc = float(has_cc.mean()) if len(df) else 0.0
    rows = []
    for fam, g in df.groupby("name", sort=True, observed=True):
        rec = family_stats(g, p, k_range)
        rec["tool"] = tool
        rec["family"] = fam
        cls = g["repeat_class_family"].mode()
        rec["class"] = str(cls.iat[0]) if len(cls) else "NA"
        rows.append(rec)
    out = pd.DataFrame(rows)
    out.attrs["frac_cons_coords"] = frac_cc
    return out


OUT_COLS = ["tool", "family", "class", "n_hits", "n_hits_conscoord", "cons_len",
            "conslen_frac", "median_len", "p90_len", "full_span_frac",
            "n_internal_endpoints", "flag", "k", "unit_len", "rho",
            "phase_enrich", "phase_z", "boundary_occ", "boundary_support",
            "best_enrich_any_k", "best_enrich_k", "reason",
            "cluster_id", "n_mates_with_len", "mate_median_cons_len",
            "L_over_mate", "unit_over_mate", "verdict"]


# --------------------------------------------------------------- cross-tool confirmation
CLUSTER_TOL = 0.5   # |L/mate_median - k| must be within this to confirm


def annotate_with_clusters(all_df: pd.DataFrame, candidates_tsv: Path,
                           tol: float = CLUSTER_TOL) -> pd.DataFrame:
    """Confirm flags against stage-0 cluster-mates -- the strongest test available.

    An over-assembly is an entry whose advertised length is ~k TIMES the length
    that independent tools give the same element.  When stage 0 has clustered
    the family with mates from other tools, that is directly checkable:

        L_over_mate = cons_len / median(cons_len of unflagged cluster-mates)

    confirmed_over_assembly  iff  |L_over_mate - k| <= tol
    Otherwise the family is suspect_coordinate_only: the block structure is
    real but the entry is not LONGER than its mates, so dividing by k would
    shorten an already-short consensus.

    Measured on GCA_951799975.1 (strict linkage): of 18 flagged families with
    length-carrying mates, exactly 4 confirm -- the four ~765 bp REPET entries
    in clusters 62, 357 and 1033, with L_over_mate 2.867 / 2.871 / 2.866 /
    2.833 against k = 3, and unit_len 254.7-257.0 vs mate lengths 264-270 bp
    (unit/mate 0.944-0.957).  The other 14 have L_over_mate 0.13-1.04: their
    consensus is SHORTER than their mates', so the block structure reflects
    fragment recruitment, not over-assembly.  This is why the pipeline should
    prefer cluster-mate inheritance over dividing by k.
    """
    df = all_df.copy()
    for c in ("cluster_id", "n_mates_with_len", "mate_median_cons_len",
              "L_over_mate", "unit_over_mate"):
        df[c] = np.nan
    df["verdict"] = ""
    df.loc[df.flag, "verdict"] = "suspect_no_cluster"

    cand = pd.read_csv(candidates_tsv, sep="\t", dtype=str)
    cand = cand[cand.cluster_id != "cluster_id"]
    member = (df.tool + ":" + df.family)
    idx = {m: i for i, m in member.items()}
    conslen = dict(zip(member, df.cons_len))
    flagged = dict(zip(member, df.flag))

    for _, row in cand.iterrows():
        mem = str(row.members).split(";")
        fl = [m for m in mem if flagged.get(m, False)]
        if not fl:
            continue
        mate_len = [conslen[m] for m in mem
                    if m not in fl and m in conslen and np.isfinite(conslen[m])]
        for m in fl:
            i = idx[m]
            df.at[i, "cluster_id"] = int(row.cluster_id)
            df.at[i, "n_mates_with_len"] = len(mate_len)
            if not mate_len:
                df.at[i, "verdict"] = "suspect_no_mate_length"
                continue
            mm = float(np.median(mate_len))
            k = float(df.at[i, "k"])
            df.at[i, "mate_median_cons_len"] = mm
            df.at[i, "L_over_mate"] = round(float(df.at[i, "cons_len"]) / mm, 3)
            df.at[i, "unit_over_mate"] = round(float(df.at[i, "unit_len"]) / mm, 3)
            df.at[i, "verdict"] = ("confirmed_over_assembly"
                                   if abs(df.at[i, "L_over_mate"] - k) <= tol
                                   else "suspect_coordinate_only")
    return df


def per_tool_table(all_df: pd.DataFrame, untestable: dict) -> pd.DataFrame:
    rows = []
    for tool, g in all_df.groupby("tool"):
        tested = g[~g.reason.isin(REJECT_REASONS)]
        fl = g[g.flag]
        rows.append(dict(
            tool=tool, status="tested",
            n_families_bed=int(len(g)), n_families_tested=int(len(tested)),
            n_flagged=int(len(fl)),
            frac_flagged=round(len(fl) / max(len(tested), 1), 4),
            modal_k=int(fl.k.mode().iat[0]) if len(fl) else 0,
            median_cons_len=float(fl.cons_len.median()) if len(fl) else float("nan"),
            median_unit_len=float(fl.unit_len.median()) if len(fl) else float("nan"),
            n_multi_boundary=int((fl.boundary_support >= 2).sum()) if len(fl) else 0))
    for tool, why in untestable.items():
        rows.append(dict(tool=tool, status=why, n_families_bed=0,
                         n_families_tested=0, n_flagged=0,
                         frac_flagged=float("nan"), modal_k=0,
                         median_cons_len=float("nan"),
                         median_unit_len=float("nan"), n_multi_boundary=0))
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--bed-dir", default="~/Documents/VGP_TEbed/inputs")
    ap.add_argument("--tools", default="rm2,edta,pantera,fastltr,repet")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache-dir", default=None,
                    help="parquet cache of the per-tool column pull")
    ap.add_argument("--eligible-tsv", default=None,
                    help="optional TSV with columns tool,family restricting the scan")
    ap.add_argument("--decoy", action="store_true",
                    help="also evaluate the half-integer decoy null and print "
                         "the specificity table")
    args = ap.parse_args(argv)

    bed_dir = Path(args.bed_dir).expanduser()
    cache_dir = Path(args.cache_dir).expanduser() if args.cache_dir else None
    elig = None
    if args.eligible_tsv:
        e = pd.read_csv(args.eligible_tsv, sep="\t")
        elig = {t: set(g.family.astype(str)) for t, g in e.groupby("tool")}

    frames, untestable, decoy_frames = [], {}, []
    for tool in args.tools.split(","):
        bed = bed_dir / f"{tool}.bed"
        if not bed.exists():
            untestable[tool] = "absent"
            continue
        cache = (cache_dir / f"{tool}_hits.parquet") if cache_dir else None
        sub = elig.get(tool) if elig else None
        df = scan_tool(tool, bed, PARAMS, cache, sub)
        if df.attrs.get("frac_cons_coords", 0.0) <= 0.0:
            untestable[tool] = "UNTESTABLE_no_cons_coords"
            print(f"[overassembly] {tool}: UNTESTABLE -- 0% of hits carry "
                  f"consensus coordinates", file=sys.stderr)
            continue
        frames.append(df)
        n_t = int((~df.reason.isin(REJECT_REASONS)).sum())
        print(f"[overassembly] {tool}: {n_t} testable families, "
              f"{int(df.flag.sum())} flagged", file=sys.stderr)
        if args.decoy:
            decoy_frames.append(scan_tool(tool, bed, PARAMS, cache, sub,
                                          k_range=DECOY_K))

    all_df = (pd.concat(frames, ignore_index=True)[OUT_COLS]
              if frames else pd.DataFrame(columns=OUT_COLS))
    out = Path(args.out).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    all_df.sort_values(["flag", "phase_enrich"], ascending=False).to_csv(
        out, sep="\t", index=False)
    print(f"[overassembly] wrote {out} ({len(all_df)} rows)", file=sys.stderr)

    tbl = per_tool_table(all_df, untestable)
    print(tbl.to_string(index=False), file=sys.stderr)

    if args.decoy and decoy_frames:
        dd = pd.concat(decoy_frames, ignore_index=True)
        n_t = int((~all_df.reason.isin(REJECT_REASONS)).sum())
        print(f"\n[overassembly] DECOY NULL (k' in {DECOY_K}) over {n_t} "
              f"testable families:", file=sys.stderr)
        print(f"  real  k in {PARAMS.k_range}: {int(all_df.flag.sum())} flagged "
              f"({all_df.flag.sum()/max(n_t,1):.4f})", file=sys.stderr)
        print(f"  decoy k' non-integer      : {int(dd.flag.sum())} flagged "
              f"({dd.flag.sum()/max(n_t,1):.4f})", file=sys.stderr)
        print(f"  implied FDR ceiling for the integer-k claim: "
              f"{dd.flag.sum()/max(int(all_df.flag.sum()),1):.2f}",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

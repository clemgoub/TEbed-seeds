#!/usr/bin/env python3
"""Sequence-free detection of OVER-ASSEMBLED library consensi.

An over-assembled consensus is one library entry that is really *k* tandem
copies of a shorter element collapsed together (RESUME.md section 5,
PIPELINE_FINDINGS.md F4).  Treating its length as the element length makes
every genuine copy look ~(1 - 1/k) truncated and sinks the family at any
near-full-length gate.

The detector needs NO library FASTA.  It uses only BED16 columns:

    implied consensus length  L = repeat_end + repeat_left   (MODAL, not max
                                  -- see stage1.consensus_length / F6)
    consensus match interval  [repeat_start, repeat_end]
    genomic hit length        chromEnd - chromStart

------------------------------------------------------------------ the rule
A family is flagged for a small integer k in K_RANGE (default 2..5) when BOTH:

(a) LENGTH COMMENSURABILITY.  The median genomic hit length m satisfies

        RHO_LO <= rho_k <= RHO_HI ,      rho_k = m * k / L

    i.e. the typical copy is one k-th of the advertised consensus.
    Defaults RHO_LO=0.70, RHO_HI=1.30.  Justification (measured):
      * upper bound -- if copies really are single units they cannot be much
        LONGER than a unit; 30% headroom absorbs indels, alignment slop and
        the fact that the modal L can under-shoot the true entry length.
      * lower bound -- the measured positives (REPET cluster 62) sit at
        rho_3 = 0.88 and 0.87; a 0.70 floor still admits families ~15% more
        degraded than the measured instance, while excluding families whose
        copies are sub-unit fragments (for those, m no longer identifies k).
      * the same two families score rho_3 = 2.50 / 2.69 for the rm2 and
        pantera consensi of the SAME element (269 / 264 bp), i.e. the window
        separates the known positive from the known negative by ~3x.
    Condition (a) on its own is weak -- every heavily-degraded family has a
    small m/L -- which is exactly why (b) carries the specificity.

(b) BLOCK STRUCTURE, not uniform spread.  Statistic = *phase-0 enrichment of
    internal match endpoints*.

      1. Endpoint pool E: every repeat_start and every repeat_end that lies
         strictly INSIDE the consensus, i.e. in (tau, L-tau) with
         tau = max(TAU_ABS, TAU_FRAC*L)  (defaults 20 bp, 5%).  Terminal
         endpoints are dropped because every family piles endpoints at 0 and
         at L -- those carry no information about internal structure and
         would score phase 0 for every k.
      2. For period P = L/k, phase phi = (x mod P)/P, binned into NB=10 bins
         with bin 0 CENTRED on phase 0 (i.e. on the block boundaries
         j*L/k, j = 1..k-1; all internal boundaries fold onto phase 0).
      3. phase_enrich_k = h[0] / (n/NB)  -- observed over expected mass in
         the boundary bin.
         phase_z_k     = (h[0] - n/NB) / sqrt(n * (1/NB) * (1 - 1/NB)).

    NULL.  Under "internal endpoints carry no k-periodic structure" the phase
    of an endpoint is uniform on [0,1), so h[0] ~ Binomial(n, 1/NB),
    E[phase_enrich] = 1 and phase_z is standard normal.  The binomial null is
    ANTI-conservative in practice (real endpoint distributions are lumpy for
    reasons unrelated to k), so the threshold is set on the *effect size*
    rather than on p:

        flag if  phase_enrich_k >= ENRICH_MIN (1.5)  AND  phase_z_k >= Z_MIN (5)

    The phase_z floor only guards small n; ENRICH_MIN does the discrimination.
    Bin 0 is anchored on phase 0 rather than taken as the maximum bin on
    purpose: an offset-free "strongest phase" version fires on any family with
    a single internal endpoint hotspot (measured: rm2 rnd-1_family-286 reaches
    2.57 that way at k=2).  Anchoring demands the hotspots sit at the
    consensus-relative positions j*L/k, which is what "the entry is k units"
    actually means.

    Measured separation on the known instance (cluster 62):
        repet:Gnig_TEdenovoGr-B-G1303-Map20  L=764  enrich_3 = 2.84  z = +39.6
        repet:Gnig_TEdenovoGr-B-G1473-Map8   L=765  enrich_3 = 2.65  z = +35.6
        rm2:rnd-1_family-286                 L=269  max over k = 0.75 (z<0)
        pantera:Unknown_572-fGobNig          L=264  max over k = 0.86 (z<0)
        rm2:rnd-5_family-3676                L=2351 max over k = 0.26 (z<0)
    No negative reaches 1.0; both positives clear 2.6.  k=3 is the unique k
    with enrichment > 1 for either positive.

If several k qualify, the reported k is the one with the largest
phase_enrich.

------------------------------------------------------------------- guards
  * MIN_HITS (50) hits carrying consensus coordinates.
  * MIN_CONSLEN (100 bp): below this, L/k is under ~20-50 bp and coordinate
    rounding dominates the phase.
  * conslen_frac >= CONSLEN_FRAC_MIN (0.5): fraction of hits whose implied
    length is within 1% of the modal L.  A family whose implied length is not
    self-consistent has no trustworthy L to divide.
  * A tool with no consensus coordinates at all (EDTA: 0/1,394,534 hits) is
    reported as UNTESTABLE, never as zero flagged.

--------------------------------------------------------------------- usage
    python tools/detect_overassembly.py \
        --bed-dir ~/Documents/VGP_TEbed/inputs \
        --tools rm2,edta,pantera,fastltr,repet \
        --out runs/GCA_951799975.1/overassembly.tsv

    # importable
    from detect_overassembly import scan_tool, family_stats, PARAMS
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, asdict
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


@dataclass(frozen=True)
class Params:
    k_range: tuple = (2, 3, 4, 5)
    n_bins: int = 10          # phase bins per period; bin 0 centred on phase 0
    tau_abs: float = 20.0     # terminal exclusion, bp
    tau_frac: float = 0.05    # terminal exclusion, fraction of L
    rho_lo: float = 0.70      # condition (a) window on median_len * k / L
    rho_hi: float = 1.30
    enrich_min: float = 1.5   # condition (b) effect size
    z_min: float = 5.0        # condition (b) significance floor
    min_hits: int = 50
    min_conslen: float = 100.0
    conslen_frac_min: float = 0.50


PARAMS = Params()


# --------------------------------------------------------------------------- core stats
def modal_consensus_length(repeat_end: np.ndarray,
                           repeat_left: np.ndarray) -> tuple[float, float]:
    """MODAL implied consensus length and the fraction of hits supporting it.

    Modal, not max: a handful of outlier rows (families flagged
    len_inconsistent upstream) imply lengths an order of magnitude too large
    -- rm2 rnd-1_family-286: modal 269 bp from 31,387 hits vs max 4,202 bp
    from 11 hits (PIPELINE_FINDINGS F6).
    """
    implied = repeat_end + repeat_left
    implied = implied[np.isfinite(implied)]
    if implied.size == 0:
        return float("nan"), float("nan")
    vals, cts = np.unique(np.round(implied).astype(np.int64), return_counts=True)
    L = float(vals[cts.argmax()])
    frac = float(np.abs(implied - L).mean() <= 0) if L == 0 else \
        float((np.abs(implied - L) <= max(1.0, 0.01 * L)).mean())
    return L, frac


def phase_enrichment(endpoints: np.ndarray, L: float, k: int,
                     n_bins: int) -> tuple[float, float, int, np.ndarray]:
    """Enrichment of internal endpoints in the phase bin centred on j*L/k.

    Returns (enrich, z, n, histogram).  Null: phi ~ U[0,1) so the boundary bin
    holds Binomial(n, 1/n_bins) of the mass.
    """
    n = int(endpoints.size)
    if n == 0:
        return float("nan"), float("nan"), 0, np.zeros(n_bins, dtype=int)
    P = L / k
    phi = ((endpoints % P) / P + 0.5 / n_bins) % 1.0      # bin 0 centred on 0
    h = np.bincount((phi * n_bins).astype(int), minlength=n_bins)[:n_bins]
    exp = n / n_bins
    enrich = float(h[0] / exp)
    z = float((h[0] - exp) / np.sqrt(exp * (1.0 - 1.0 / n_bins)))
    return enrich, z, n, h


def boundary_occupancy(endpoints: np.ndarray, L: float, k: int,
                       half_width_frac: float = 0.05) -> list[float]:
    """Share of internal endpoints within +-half_width_frac*P of each j*L/k."""
    P = L / k
    w = half_width_frac * P
    n = max(endpoints.size, 1)
    return [round(float(np.abs(endpoints - j * P).__le__(w).sum()) / n, 3)
            for j in range(1, k)]


def family_stats(hits: pd.DataFrame, p: Params = PARAMS) -> dict:
    """Full per-family record.  `hits` needs repeat_start/end/left + hit_len."""
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
               n_internal_endpoints=0, k=0, rho=float("nan"),
               phase_enrich=float("nan"), phase_z=float("nan"),
               best_enrich_any_k=float("nan"), best_enrich_k=0,
               boundary_occ="", flag=False, reason="")

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

    tau = max(p.tau_abs, p.tau_frac * L)
    pool = np.concatenate([rs, re_])
    E = pool[(pool > tau) & (pool < L - tau)]
    rec["n_internal_endpoints"] = int(E.size)
    if E.size < p.min_hits:
        rec["reason"] = "too_few_internal_endpoints"
        return rec

    m = float(np.median(hl))
    best_any = (-np.inf, 0)
    passing = []
    for k in p.k_range:
        enrich, z, _, _ = phase_enrichment(E, L, k, p.n_bins)
        if enrich > best_any[0]:
            best_any = (enrich, k)
        rho = m * k / L
        cond_a = p.rho_lo <= rho <= p.rho_hi
        cond_b = (enrich >= p.enrich_min) and (z >= p.z_min)
        if cond_a and cond_b:
            passing.append((enrich, k, rho, z))
    rec["best_enrich_any_k"] = round(float(best_any[0]), 3)
    rec["best_enrich_k"] = int(best_any[1])

    if not passing:
        rec["reason"] = "no_k_passes"
        return rec

    enrich, k, rho, z = max(passing)
    rec.update(k=int(k), rho=round(float(rho), 3),
               phase_enrich=round(float(enrich), 3), phase_z=round(float(z), 2),
               boundary_occ=";".join(str(v) for v in boundary_occupancy(E, L, k)),
               flag=True, reason="over_assembly")
    return rec


# --------------------------------------------------------------------------- driver
def load_tool(bed: Path, cache: Path | None = None) -> pd.DataFrame:
    """Read only the BED16 columns the detector needs (the BEDs are 120-240 MB)."""
    if cache is not None and cache.exists():
        return pd.read_parquet(cache)
    parts = []
    for ch in pd.read_csv(bed, sep="\t", names=BED16_COLS, skiprows=1,
                          usecols=USECOLS, na_values=["NA"], chunksize=2_000_000,
                          dtype={"name": "string", "repeat_class_family": "string"}):
        for c in ("repeat_start", "repeat_end", "repeat_left"):
            ch[c] = pd.to_numeric(ch[c], errors="coerce")
        ch["hit_len"] = ch["chromEnd"] - ch["chromStart"]
        parts.append(ch[["name", "repeat_class_family", "repeat_start",
                         "repeat_end", "repeat_left", "hit_len"]])
    df = pd.concat(parts, ignore_index=True)
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(cache, index=False)
    return df


def scan_tool(tool: str, bed: Path, p: Params = PARAMS,
              cache: Path | None = None,
              eligible: set[str] | None = None) -> pd.DataFrame:
    """Per-family records for one tool BED.  Empty frame if no cons coords."""
    df = load_tool(bed, cache)
    df["name"] = df["name"].astype(str)
    if eligible is not None:
        df = df[df["name"].isin(eligible)]
    has_cc = df["repeat_end"].notna() & df["repeat_left"].notna()
    frac_cc = float(has_cc.mean()) if len(df) else 0.0
    rows = []
    for fam, g in df.groupby("name", sort=True, observed=True):
        rec = family_stats(g, p)
        rec["tool"] = tool
        rec["family"] = fam
        cls = g["repeat_class_family"].mode()
        rec["class"] = str(cls.iat[0]) if len(cls) else "NA"
        rows.append(rec)
    out = pd.DataFrame(rows)
    out.attrs["frac_cons_coords"] = frac_cc
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bed-dir", default="~/Documents/VGP_TEbed/inputs")
    ap.add_argument("--tools", default="rm2,edta,pantera,fastltr,repet")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache-dir", default=None,
                    help="parquet cache of the per-tool column pull")
    ap.add_argument("--eligible-tsv", default=None,
                    help="optional 2-col TSV (tool, family) restricting the scan")
    args = ap.parse_args(argv)

    bed_dir = Path(args.bed_dir).expanduser()
    cache_dir = Path(args.cache_dir).expanduser() if args.cache_dir else None
    elig = None
    if args.eligible_tsv:
        e = pd.read_csv(args.eligible_tsv, sep="\t")
        elig = {t: set(g.family) for t, g in e.groupby("tool")}

    frames, summary = [], []
    for tool in args.tools.split(","):
        bed = bed_dir / f"{tool}.bed"
        if not bed.exists():
            summary.append(dict(tool=tool, status="absent", n_families_tested=0,
                                n_flagged=0))
            continue
        cache = (cache_dir / f"{tool}_hits.parquet") if cache_dir else None
        df = scan_tool(tool, bed, PARAMS, cache,
                       elig.get(tool) if elig else None)
        fcc = df.attrs.get("frac_cons_coords", 0.0)
        if fcc <= 0.0:
            summary.append(dict(tool=tool, status="UNTESTABLE_no_cons_coords",
                                n_families_tested=0, n_flagged=0))
            print(f"[overassembly] {tool}: UNTESTABLE -- 0% of hits carry "
                  f"consensus coordinates", file=sys.stderr)
            continue
        frames.append(df)
        tested = df[df.reason != "too_few_hits_with_cons_coords"]
        tested = tested[~tested.reason.isin(["consensus_too_short",
                                             "consensus_length_inconsistent",
                                             "too_few_internal_endpoints"])]
        summary.append(dict(tool=tool, status="tested",
                            n_families_tested=int(len(tested)),
                            n_flagged=int(df.flag.sum())))
        print(f"[overassembly] {tool}: {len(tested)} testable families, "
              f"{int(df.flag.sum())} flagged", file=sys.stderr)

    out = Path(args.out).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["tool", "family", "class", "n_hits", "n_hits_conscoord", "cons_len",
            "conslen_frac", "median_len", "p90_len", "n_internal_endpoints",
            "flag", "k", "rho", "phase_enrich", "phase_z", "boundary_occ",
            "best_enrich_any_k", "best_enrich_k", "reason"]
    all_df = (pd.concat(frames, ignore_index=True)[cols]
              if frames else pd.DataFrame(columns=cols))
    all_df["unit_len"] = np.where(all_df.flag, all_df.cons_len / all_df.k.replace(0, np.nan),
                                  np.nan).round(1)
    all_df.sort_values(["flag", "phase_enrich"], ascending=False).to_csv(
        out, sep="\t", index=False)
    print(f"[overassembly] wrote {out} ({len(all_df)} rows)", file=sys.stderr)
    print(pd.DataFrame(summary).to_string(index=False), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Protein-domain evidence for a rebuilt consensus, and the LHF size/domain gate.

Why this exists: the queue was ranked by `pooled_full_len`, which is a COPY
NUMBER ranking, so it promoted abundant short non-autonomous elements. The
top-20 had a median rebuilt consensus of 469 bp with 19 of 20 under 1 kb --
MITE and SINE scale, not low-hanging fruit for a curated library.

The evidence is cheap enough to run on every packet: 235 consensi (0.50 Mb)
scan in 1.8 s with BATH against TE-Aid's 130 curated TE-protein Pfam profiles.
There is no reason to treat it as a one-off.

TWO CORRECTIONS TO THE FIRST DRAFT OF THIS GATE, both of which changed the
answer rather than refining it:

1. The core-transposase list omitted **DDE_3** (PF13358, the Tc1/Mariner
   catalytic domain) and **HTH_Tnp_Tc3_1/2** (its DNA-binding domain). Tc1/
   Mariner is among the most abundant DNA transposon superfamilies in teleost
   fish, and a genome-wide scan of this assembly finds 1,629 transposase-domain
   loci dominated by exactly those two (HTH_Tnp_Tc3_2 682, DDE_3 401). The
   original gate returned ZERO tier-A DNA transposons, and that zero was this
   omission, not the genome.
2. `Transposase_22` was on the transposase side but is the LINE L1 ORF-
   associated profile; it made an rm2 `#LINE/L1-Tx1` family read as a DNA
   transposon. It belongs with the accessory profiles.

The size floors carry the same hazard and are deliberately generous. Autonomous
Tc1/Mariner in this genome measures 1,076-1,633 bp, so a flat DNA floor of
1,500 bp rejects real autonomous elements as "fragment scale". A floor that
manufactures its own conclusion is worse than no floor.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------- profile sets
# Core replication/mobilisation machinery: possessing one of these is what makes
# an element autonomous. Grouped by the enzyme, because "which core domain"
# is more informative to a curator than a single boolean.
CORE_RT = (
    r"RVT_\d+$", r"RVT_N$", r"Pao_retrotransp$", r"RT_RNaseH", r"Exo_endo_phos",
    r"TERT$", r"RVP", r"Peptidase_A17$",
)
CORE_TRANSPOSASE = (
    # DDE catalytic cores
    r"DDE_Tnp_", r"DDE_\d+$", r"Transposase_(?!22$)\d+$", r"Dimer_Tnp_hAT$",
    r"DBD_Tnp_", r"MULE$", r"Mutator", r"Plant_tran$", r"Tnp_P_element",
    # Tc1/Mariner: the catalytic DDE_3 and its HTH DNA-binding domains.
    # Omitting these is what produced a queue with no DNA transposons.
    r"HTH_Tnp_Tc3_", r"HTH_Tnp_Tc5$", r"Transp_Tc5_C$", r"HTH_Tnp_4$",
    r"HTH_48$", r"zf-Tnp_2$",
    # rolling-circle (Helitron) and Replitron replication initiators
    r"Helitron_like_N$", r"Rol_Rep_N$", r"Replitron_HUH$", r"Rep_3$",
)
CORE_INTEGRASE = (r"^rve", r"Integrase_H2C2$", r"IN_DBD_C$")

CORE_SETS = {"RT": CORE_RT, "transposase": CORE_TRANSPOSASE,
             "integrase": CORE_INTEGRASE}

# Per-order minimum consensus length for an element to be plausibly autonomous.
# SINEs are exempt (non-autonomous by definition) but are capped instead: an
# oversized SINE is a chimera or an rDNA array, not a SINE.
MIN_LEN = {"LTR": 4000, "DIRS": 4000, "LINE": 2500, "RC": 3000,
           "DNA": 1100, "TIR": 1100, "PLE": 1200, "Helitron": 3000}
SINE_MAX = 700

MIN_HMM_COV = 0.50      # fraction of the profile the hit must cover
MIN_ALN_ROWS = 10       # tier-A precondition: a 5-row seed is not a consensus
MIN_DEPTH = 5


def _matches(name: str, patterns) -> bool:
    return any(re.search(p, name) for p in patterns)


def classify_domain(name: str) -> str | None:
    """'RT' | 'transposase' | 'integrase' for a core domain, else None."""
    for kind, pats in CORE_SETS.items():
        if _matches(name, pats):
            return kind
    return None


# ------------------------------------------------------------------ scanning
def run_bath(fasta: Path, out_tbl: Path, bath_bin: Path, profiles: Path,
             cpu: int = 8, evalue: float = 1e-5) -> dict:
    """BATH: translated search of protein profiles against DNA sequences.

    BATH is frameshift-aware, which matters here -- the consensi are rebuilt
    from degraded genomic copies and a plain six-frame translation loses
    domains to indels.
    """
    out_tbl.parent.mkdir(parents=True, exist_ok=True)
    cmd = [str(bath_bin), "--cpu", str(cpu), "-E", str(evalue),
           "--tblout", str(out_tbl), str(profiles), str(fasta)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    return dict(cmd=" ".join(cmd), returncode=p.returncode,
                stderr_tail=(p.stderr or "")[-1500:])


def parse_bath_tbl(path: Path) -> pd.DataFrame:
    """BATH --tblout -> tidy hits with HMM coverage.

    The BATH table carries a leading hit-ID column, so the fields are shifted
    one right of the HMMER layout they otherwise resemble:

        1  target-name  acc  query-name  acc  hmm_len  hmm_from  hmm_to
           seq_len  ali_from  ali_to  E-value  score  bias  PID  desc
        0       1       2         3      4       5        6        7
                  8         9        10      11      12    13   14

    Reading them one column left silently produced ZERO domain hits from a
    495-hit table, which then read as "no autonomous elements in this genome".
    """
    cols = ["seq", "profile", "hmm_len", "hmm_from", "hmm_to", "evalue"]
    rows = []
    if not Path(path).exists():
        return pd.DataFrame(columns=cols + ["hmm_cov"])
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 12:
            continue
        try:
            seq, profile = f[1], f[3]
            hmm_len, hf, ht = int(f[5]), int(f[6]), int(f[7])
            ev = float(f[11])
        except (ValueError, IndexError):
            continue
        hmm_cov = (abs(ht - hf) + 1) / hmm_len if hmm_len else float("nan")
        rows.append(dict(seq=seq, profile=profile, hmm_len=hmm_len,
                         hmm_from=hf, hmm_to=ht, evalue=ev, hmm_cov=hmm_cov))
    return pd.DataFrame(rows)


def summarise_domains(hits: pd.DataFrame, min_cov: float = MIN_HMM_COV) -> pd.DataFrame:
    """Per sequence: which core domains it carries, at what coverage."""
    if not len(hits):
        return pd.DataFrame(columns=["seq", "n_dom", "domains", "core_kinds",
                                     "n_core", "best_cov", "core_domains"])
    h = hits.copy()
    h["kind"] = h.profile.map(classify_domain)
    strong = h[h["hmm_cov"] >= min_cov]   # NOT h.cov: that is DataFrame.cov()
    out = []
    for seq, g in h.groupby("seq"):
        sg = strong[strong.seq == seq]
        core = sg[sg.kind.notna()]
        out.append(dict(
            seq=seq, n_dom=int(g.profile.nunique()),
            domains=";".join(sorted(g.profile.unique())[:12]),
            core_kinds=";".join(sorted(core.kind.dropna().unique())),
            n_core=int(core.profile.nunique()),
            core_domains=";".join(sorted(core.profile.unique())),
            best_cov=float(g["hmm_cov"].max()) if len(g) else float("nan"),
        ))
    return pd.DataFrame(out)


# ---------------------------------------------------------------------- gate
def order_of(canonical_path, tp) -> str:
    """Coarse order for the size rule, from the canonical path or the TP."""
    s = f"{canonical_path or ''}|{tp or ''}"
    for key, pat in (("SINE", r"SINE"), ("LINE", r"LINE"), ("PLE", r"PLE|Penelope"),
                     ("DIRS", r"DIRS"), ("LTR", r"LTR"),
                     ("RC", r"Helitron|RC\b|Rolling"), ("TIR", r"TIR|DNA")):
        if re.search(pat, s, re.I):
            return key
    return "Unknown"


def assign_tier(row) -> str:
    """A-F, on consensus length, core-domain evidence and alignment support."""
    order, bp = row["order"], row["bp"]
    has_core = row.get("n_core", 0) > 0
    has_any = row.get("n_dom", 0) > 0
    floor = MIN_LEN.get(order)
    thin = (row.get("aln_rows", 0) < MIN_ALN_ROWS
            or row.get("median_depth", 0) < MIN_DEPTH)

    if order == "SINE":
        if bp > SINE_MAX:
            return "D - oversized for class, no core domain (suspect chimera/rDNA)"
        return "E - non-autonomous at expected size (SINE)"
    if has_core and floor and bp >= floor:
        return ("B - core domain, alignment too thin to trust" if thin
                else "A - autonomous (core domain + class size)")
    if has_core:
        return "B - core domain, short for class"
    if has_any:
        return "C - accessory domain only (likely truncated)"
    if floor and bp >= floor:
        return "D - oversized for class, no core domain (suspect chimera/rDNA)"
    return "F - too short for class (MITE/fragment scale)"


def rank_tier_a(df: pd.DataFrame) -> pd.DataFrame:
    """Order tier A by distinct core domains, then near-full-length copies.

    `pooled_full_len` alone is a copy-number ranking and promotes abundant
    short elements; a complete Gypsy carrying six core domains at full HMM
    coverage is the better head of an LHF queue.
    """
    a = df[df.tier.str.startswith("A")].copy()
    return a.sort_values(["n_core", "best_cov", "n_full", "bp"],
                         ascending=[False, False, False, False])

"""Dfam-compliant Stockholm seed emission.

Format authority: vendor/dfam-curator/Dfam_Seeds.md, validated with `stk lint`.
The rules that actually bite, all of them checked by lint:

  - DE, AU, TP, OC and SQ are REQUIRED, plus a `#=GC RF` line
    (`missing_required_field` / `rf_missing`, both ERROR).
  - Dfam standardises on `.` as the gap character; `-` is read but warned
    (`seq_nonstandard_gap`, one warning per row).
  - One line per sequence.  Block/interleaved format is an ERROR
    (`block_format`), so no wrapping however wide the alignment gets.
  - SQ must equal the actual row count (`sq_mismatch`).
  - DE is capped at 80 characters (`de_too_long`).
  - AU must be a real full name.  `au_format` is an ERROR for
    "Last Initial" style or a single-letter first name -- so the placeholder
    "C. Goubert" in config would fail; a spelled-out first name is required.
  - RF carries the consensus at match columns and `.` at insert columns, and
    lint compares it against its own consensus caller (`rf_consensus_mismatch`,
    WARN).  `stk edit --update-consensus` rewrites RF with Dfam's own caller,
    which is the reliable way to agree with it.

Identifiers are Smitten: assembly:sequence:start-end_strand, 1-based fully
closed.  Public accessions matter -- Dfam's TSD and extension algorithms
retrieve flanking sequence by them, so an opaque local label makes a seed
unusable downstream even if it lints clean.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np

GAP = "."
REQUIRED = ("DE", "AU", "TP", "OC", "SQ")


def rf_line(cons_chars: str, is_match: np.ndarray) -> str:
    """Full-width RF: consensus base at match columns, '.' at insert columns."""
    out = []
    it = iter(cons_chars)
    for m in is_match:
        out.append(next(it) if m else GAP)
    return "".join(out)


def format_record(ids: list[str], rows: list[str], meta: dict) -> str:
    """One Stockholm record.  `meta` keys are #=GF tags; list values repeat the
    tag (OC and CC legitimately occur several times)."""
    if len(ids) != len(rows):
        raise ValueError("ids/rows length mismatch")
    widths = {len(r) for r in rows}
    if len(widths) > 1:
        raise ValueError(f"ragged alignment: widths {sorted(widths)}")
    if "RF" in meta and len(meta["RF"]) != rows[0].__len__():
        raise ValueError("RF width does not match the alignment")

    lines = ["# STOCKHOLM 1.0"]
    order = ["ID", "DE", "AU", "SE", "TP", "OC", "SQ", "TD", "CT", "KD", "BM",
             "RN", "RM", "RD", "DR", "CC", "**"]
    for tag in order:
        if tag not in meta or tag == "RF":
            continue
        val = meta[tag]
        for v in (val if isinstance(val, (list, tuple)) else [val]):
            lines.append(f"#=GF {tag:<5s} {v}")
    pad = max(len(i) for i in ids) + 4
    if "RF" in meta:
        lines.append(f"{'#=GC RF':<{pad}s}{meta['RF']}")
    for i, r in zip(ids, rows):
        lines.append(f"{i:<{pad}s}{r}")
    if "MM" in meta:
        lines.append(f"{'#=GC MM':<{pad}s}{meta['MM']}")
    lines.append("//")
    return "\n".join(lines) + "\n"


def write_stockholm(path: Path, records: list[str]) -> None:
    with open(path, "w") as fh:
        for r in records:
            fh.write(r)


def to_rows(m: np.ndarray, lo: int, hi: int) -> list[str]:
    """Alignment matrix -> Stockholm rows, with Dfam's '.' gap character."""
    return ["".join(chr(c) for c in row).replace("-", GAP)
            for row in m[:, lo:hi]]


def lint(stk: Path, stk_bin: Path, genome: Path | None = None,
         no_network: bool = True, min_severity: str = "info") -> dict:
    cmd = [str(stk_bin), "lint", "--min-severity", min_severity]
    if no_network:
        cmd.append("--no-network")
    if genome:
        cmd += ["--genome", str(genome)]
    cmd.append(str(stk))
    p = subprocess.run(cmd, capture_output=True, text=True)
    out = (p.stdout or "") + (p.stderr or "")
    codes: dict[str, int] = {}
    for line in out.splitlines():
        for sev in ("ERROR", "WARN", "INFO"):
            if sev in line:
                # lines look like "... ERROR  code_name  detail"
                parts = line.split()
                if sev in parts:
                    j = parts.index(sev)
                    if j + 1 < len(parts):
                        codes[f"{sev}:{parts[j + 1]}"] = \
                            codes.get(f"{sev}:{parts[j + 1]}", 0) + 1
                break
    return dict(cmd=" ".join(cmd), returncode=p.returncode, output=out,
                codes=codes,
                n_error=sum(v for k, v in codes.items() if k.startswith("ERROR")))


def update_consensus(stk: Path, out: Path, stk_bin: Path) -> bool:
    """Rewrite #=GC RF with Dfam's own consensus caller, so lint's
    rf_consensus_mismatch check compares like with like."""
    p = subprocess.run([str(stk_bin), "edit", "--update-consensus", str(stk)],
                       capture_output=True, text=True)
    if p.returncode != 0:
        return False
    out.write_text(p.stdout)
    return True

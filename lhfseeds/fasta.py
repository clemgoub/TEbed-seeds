"""Random-access FASTA reader over a samtools .fai index.

The assembly is ~880 MB; loading it into memory per cluster is wasteful and
Bio.SeqIO.index re-parses the whole file. The .fai already carries everything
needed for an O(1) seek, so use it directly.
"""
from __future__ import annotations

from pathlib import Path

_COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn",
                      "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


class IndexedFasta:
    """Minimal faidx-backed reader. fetch() is 0-based half-open, like BED."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        fai = self.path.with_suffix(self.path.suffix + ".fai")
        if not fai.exists():
            raise FileNotFoundError(f"missing index {fai}; run `samtools faidx {path}`")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        for line in open(fai):
            name, length, offset, line_bases, line_width = line.split("\t")[:5]
            self.index[name] = (int(length), int(offset), int(line_bases),
                                int(line_width))
        self._fh = open(self.path, "rb")

    def __contains__(self, chrom: str) -> bool:
        return chrom in self.index

    def length(self, chrom: str) -> int:
        return self.index[chrom][0]

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """0-based half-open [start, end). Clamped to the contig."""
        length, offset, lb, lw = self.index[chrom]
        start = max(0, start)
        end = min(length, end)
        if end <= start:
            return ""
        # byte offset of a base = offset + (pos // line_bases) * line_width
        #                                + (pos % line_bases)
        b0 = offset + (start // lb) * lw + (start % lb)
        b1 = offset + ((end - 1) // lb) * lw + ((end - 1) % lb) + 1
        self._fh.seek(b0)
        raw = self._fh.read(b1 - b0).decode("ascii")
        return "".join(raw.split())

    def close(self) -> None:
        self._fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

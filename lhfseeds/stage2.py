"""Stage 2 -- from copy coordinates to extracted sequence and an MSA.

Three things happen here, in order:

1. LOCUS DEDUPLICATION.  Stage 1 emits one row per (member, copy).  Members of
   a cluster are, by construction, different tools describing the SAME element,
   so they annotate the same genomic loci over and over: measured on cluster 62,
   90,313 pooled copy rows collapse to 34,190 distinct loci -- 2.64 annotations
   per locus, up to 30.  Sampling the pooled rows naively would put the same
   sequence into the alignment ~2.6x on average, inflating apparent depth and
   biasing the consensus toward whichever tool is most prolific.  One
   representative per locus, with the cross-tool support recorded.

2. FRAGMENT EXPLOSION.  A copy may be several fragments bridged across a gap.
   The gap interior is excluded on purpose (62.9% of rm2 inter-fragment gaps
   are mostly another element).  Emitting one alignment row per FRAGMENT keeps
   every Dfam sequence identifier a genuine contiguous interval while still
   excluding the gap -- both requirements at once, no concatenated chimeras.

3. EXTRACTION.  Two files per packet, because they serve different consumers:
     copies.fa          element only          -> MSA, seed, consensus
     copies.flanked.fa  element +- flank_bp   -> TSD detection, TE-Aid, Refiner
   Seed identifiers therefore describe the element interval exactly, which is
   what `stk lint --genome` validates against the assembly.
"""
from __future__ import annotations

import collections
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from .fasta import IndexedFasta, revcomp

LOCUS_MIN_RECIPROCAL = 0.5   # same rule stage0 uses for graph edges
MIN_FRAG_BP = 30             # below this a fragment carries no alignable signal
MIN_MATCH_DEPTH = 3          # a match column needs this many rows reaching it


# --------------------------------------------------------------- 1. loci
def _locus_ids(starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """Group intervals on one chromosome into loci.

    Two annotations are the same locus when they overlap reciprocally by
    >= LOCUS_MIN_RECIPROCAL of the SHORTER one.  Mere adjacency is not enough:
    this element occurs in tandem arrays (2,892 same-strand runs measured), and
    a touch-based rule would chain a whole array into a single "locus" -- the
    same failure mode as F1 one level up.
    """
    n = len(starts)
    parent = np.arange(n)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return int(x)

    order = np.argsort(starts, kind="stable")
    frontier: list[int] = []          # intervals still able to overlap
    for oi in order:
        s, e = starts[oi], ends[oi]
        keep = []
        for mi in frontier:
            ms, me = starts[mi], ends[mi]
            if me <= s:
                continue              # cannot overlap this or anything later
            keep.append(mi)
            inter = min(e, me) - max(s, ms)
            if inter >= LOCUS_MIN_RECIPROCAL * min(e - s, me - ms):
                ra, rb = find(oi), find(mi)
                if ra != rb:          # union, so a chain gets ONE id, not the last
                    parent[ra] = rb
        keep.append(oi)
        frontier = keep
    roots = np.array([find(i) for i in range(n)])
    _, lid = np.unique(roots, return_inverse=True)
    return lid.astype(np.int64)


def cluster_loci(copies: pd.DataFrame, modal_cons_len: float = float("nan")
                 ) -> pd.DataFrame:
    """Collapse pooled per-member copies to one representative per locus.

    Representative ranking, in order:
      1. near-full-length             -- a complete copy beats a stub
      2. own consensus coordinates    -- bed16 beats an inherited length
      3. fewest fragments             -- an unbridged copy beats a bridged one
      4. longest span                 -- the runaway cap already bounds this

    `is_full` is pooled with OR: if any tool at the locus reached both consensus
    ends, the locus is full-length.  This is where pooling earns its keep --
    EDTA carries no consensus coordinates at all and can never call a copy full
    on its own, but at a locus rm2 also annotates it inherits that judgement.
    """
    if not len(copies):
        return copies.assign(locus_id=[], n_members=[], n_tools=[],
                             members=[], any_full=[], span_disagreement_bp=[])
    df = copies.reset_index(drop=True).copy()
    df["_tool"] = df.member.str.split(":").str[0]
    lid = np.full(len(df), -1, dtype=np.int64)
    base = 0
    for ch, g in df.groupby("chrom", sort=False):
        sub = _locus_ids(g.start.to_numpy(), g.end.to_numpy())
        lid[g.index.to_numpy()] = sub + base
        base += sub.max() + 1 if len(sub) else 0
    df["locus_id"] = lid

    df["_rank_full"] = (~df.is_full.astype(bool)).astype(int)          # 0 best
    df["_rank_src"] = (df.cons_len_source != "bed16").astype(int)
    df["_rank_frag"] = df.n_fragments.astype(int)
    df["_rank_span"] = -(df.end - df.start)

    agg = df.groupby("locus_id").agg(
        n_members=("member", "size"),
        n_tools=("_tool", "nunique"),
        members=("member", lambda s: ";".join(sorted(set(s)))),
        any_full=("is_full", "any"),
        locus_start=("start", "min"),
        locus_end=("end", "max"),
    )
    # drop_duplicates, NOT groupby().first(): first() is per-column and takes the
    # first NON-NULL value, which would splice fields from several annotations
    # into a representative that does not exist at the locus.
    rep = (df.sort_values(["locus_id", "_rank_full", "_rank_src",
                           "_rank_frag", "_rank_span"], kind="stable")
             .drop_duplicates(subset="locus_id", keep="first"))
    out = rep.merge(agg, on="locus_id")
    out["span_disagreement_bp"] = ((out.locus_end - out.locus_start)
                                   - (out.end - out.start))
    if not np.isnan(modal_cons_len):
        # a locus far wider than the element is a chained array, not one copy
        out["locus_x_consensus"] = ((out.locus_end - out.locus_start)
                                    / modal_cons_len).round(2)
    return out.drop(columns=[c for c in out.columns if c.startswith("_rank")]
                    + ["_tool"])


# ---------------------------------------------------- 2. fragment explosion
def explode_fragments(copies: pd.DataFrame, min_frag_bp: int = MIN_FRAG_BP
                      ) -> pd.DataFrame:
    """One row per fragment, each a genuine contiguous interval.

    A single-fragment copy yields exactly one row, so this is a no-op for the
    99.75% of gap-aware copies that are unbridged.
    """
    rows = []
    for r in copies.itertuples():
        fs = [int(x) for x in str(r.frag_starts).split(";")]
        fe = [int(x) for x in str(r.frag_ends).split(";")]
        nf = len(fs)
        for i, (s, e) in enumerate(zip(fs, fe)):
            if e - s < min_frag_bp:
                continue
            d = {c: getattr(r, c) for c in copies.columns
                 if c not in ("start", "end", "frag_starts", "frag_ends")}
            d.update(start=s, end=e, frag_index=i, frag_of=nf,
                     is_fragment_of_bridged=nf > 1,
                     strand=norm_strand(d.get("strand")),
                     strand_called=d.get("strand") in ("+", "-"))
            rows.append(d)
    return pd.DataFrame(rows)


# --------------------------------------------------------- 3. extraction
def norm_strand(strand) -> str:
    """BED strand -> '+' or '-'.

    EDTA leaves the strand uncalled ('.') on ~0.3% of its hits.  Extraction
    reverse-complements only on '-', so an uncalled strand IS extracted forward
    and the identifier must say '+'.  Leaving '.' in play let extraction and
    identifier retagging disagree about the orientation of the same row; the
    Smitten format has no '.' anyway, and `stk lint --genome` caught it.
    """
    return "-" if strand == "-" else "+"


def seq_id(assembly: str, chrom: str, start: int, end: int, strand: str) -> str:
    """Dfam Smitten identifier: 1-based fully closed, from BED half-open."""
    return f"{assembly}:{chrom}:{start + 1}-{end}_{strand}"


def extract(fa: IndexedFasta, rows: pd.DataFrame, assembly: str,
            flank_bp: int = 0) -> tuple[list[tuple[str, str]], pd.DataFrame]:
    """Extract each row; reverse-complement '-' rows.

    Returns (records, stats).  With flank_bp>0 the flank is added in GENOMIC
    space before the reverse complement, so a '-' copy gets its upstream flank
    where the element's 5' end actually is.
    """
    recs, stats = [], []
    for r in rows.itertuples():
        a, b = int(r.start) - flank_bp, int(r.end) + flank_bp
        s = fa.fetch(r.chrom, a, b)
        got_left = int(r.start) - max(0, a)
        got_right = min(fa.length(r.chrom), b) - int(r.end)
        if r.strand == "-":
            s = revcomp(s)
            got_left, got_right = got_right, got_left
        sid = seq_id(assembly, r.chrom, int(r.start), int(r.end), r.strand)
        recs.append((sid, s))
        n_count = s.upper().count("N")
        stats.append(dict(seq_id=sid, length=len(s), n_bases=n_count,
                          n_frac=n_count / max(len(s), 1),
                          flank_left=got_left, flank_right=got_right,
                          element_bp=int(r.end) - int(r.start)))
    return recs, pd.DataFrame(stats)


def write_fasta(recs: list[tuple[str, str]], path: Path, width: int = 60) -> None:
    with open(path, "w") as fh:
        for sid, s in recs:
            fh.write(f">{sid}\n")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    recs, sid, buf = [], None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if sid is not None:
                recs.append((sid, "".join(buf)))
            sid, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if sid is not None:
        recs.append((sid, "".join(buf)))
    return recs


# ------------------------------------------------------------ 4. sampling
def stratified_sample(copies: pd.DataFrame, cap: int, floor: int,
                      full_frac: float, seed: int = 42,
                      full_col: str = "any_full") -> pd.DataFrame:
    """Sample <=cap copies, ~full_frac of them near-full-length.

    Partial copies are stratified by WHICH PART of the consensus they cover, in
    five bands over the consensus midpoint.  Uniform random sampling of partials
    would over-represent whichever end of the element survives most often (for a
    LINE, the 3' end), and the rebuilt consensus would then be supported at one
    end and thin at the other.  Realized -- not target -- composition is what
    the packet records, because repeatome intactness varies by taxon.
    """
    rng = np.random.default_rng(seed)
    isfull = copies[full_col].astype(bool)
    if len(copies) <= max(floor, 1):
        return copies.assign(sampled_as=np.where(isfull, "full", "partial"),
                             cov_band=-1)
    n = min(cap, len(copies))
    full = copies[isfull]
    part = copies[~isfull]

    want_full = min(int(round(n * full_frac)), len(full))
    take_full = (full.sample(want_full, random_state=int(rng.integers(1 << 31)))
                 if want_full else full.iloc[:0])
    take_full = take_full.assign(cov_band=-1)

    want_part = n - len(take_full)
    if want_part > 0 and len(part):
        cs, ce, cl = part.cons_start, part.cons_end, part.cons_len
        mid = (cs + ce) / 2 / cl.replace(0, np.nan)
        band = pd.cut(mid, bins=[-np.inf, .2, .4, .6, .8, np.inf],
                      labels=False).fillna(-1).astype(int)
        part = part.assign(cov_band=band)
        picks, per = [], max(1, want_part // max(band.nunique(), 1))
        for b, g in part.groupby("cov_band"):
            picks.append(g.sample(min(per, len(g)),
                                  random_state=int(rng.integers(1 << 31))))
        take_part = pd.concat(picks) if picks else part.iloc[:0]
        if len(take_part) < want_part:            # top up from what is left
            left = part.drop(take_part.index)
            if len(left):
                take_part = pd.concat([take_part, left.sample(
                    min(want_part - len(take_part), len(left)),
                    random_state=int(rng.integers(1 << 31)))])
        take_part = take_part.iloc[:want_part]
    else:
        take_part = part.iloc[:0].assign(cov_band=-1)

    short = n - len(take_full) - len(take_part)   # not enough partials existed
    if short > 0 and len(full) > len(take_full):
        extra = full.drop(take_full.index).sample(
            min(short, len(full) - len(take_full)),
            random_state=int(rng.integers(1 << 31))).assign(cov_band=-1)
        take_full = pd.concat([take_full, extra])

    out = pd.concat([take_full.assign(sampled_as="full"),
                     take_part.assign(sampled_as="partial")])
    return out.sort_values(["chrom", "start"]).reset_index(drop=True)


# ------------------------------------------------------------------ 5. MSA
def mafft(in_fa: Path, out_fa: Path, threads: int = 0,
          extra: tuple[str, ...] = ("--auto", "--quiet",
                                    "--adjustdirectionaccurately")) -> dict:
    """Run MAFFT.  --adjustdirectionaccurately because tool strand calls
    disagree at some loci; letting MAFFT settle it costs little and a
    strand-flipped row would otherwise poison the consensus."""
    cmd = ["mafft", *extra, "--thread", str(threads or 1), str(in_fa)]
    with open(out_fa, "w") as fh:
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"mafft failed ({p.returncode}): {p.stderr[-2000:]}")
    return dict(cmd=" ".join(cmd), stderr_tail=p.stderr[-500:])


_GAP_CHARS = str.maketrans({".": "-", "_": "-", "~": "-"})


def msa_matrix(recs: list[tuple[str, str]]) -> tuple[list[str], np.ndarray]:
    """(ids, uint8 matrix), gap characters normalised to '-'.

    Stockholm permits four gap characters and different tools pick different
    ones: MAFFT writes '-', Refiner and Dfam write '.'.  Treating only '-' as a
    gap made every Refiner column look fully occupied, so trimming did nothing,
    every column became a match column, and the "consensus" came out as the
    full alignment width (349 bp against Refiner's own 267).
    MAFFT also lowercases rows it reversed and prefixes '_R_' on their ids.
    """
    ids = [i[3:] if i.startswith("_R_") else i for i, _ in recs]
    w = max(len(s) for _, s in recs)
    m = np.full((len(recs), w), ord("-"), dtype=np.uint8)
    for i, (_, s) in enumerate(recs):
        a = np.frombuffer(s.upper().translate(_GAP_CHARS).encode(), dtype=np.uint8)
        m[i, :len(a)] = a
    return ids, m


def span_occupancy(m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-column (occupancy, spanning-row count), normalised by ROWS THAT
    REACH THE COLUMN rather than by every row in the alignment.

    A copy that stops at position 900 says nothing about position 3000, so
    counting it in position 3000's denominator penalises a column for evidence
    that was never available. Measured cost of getting this wrong: on the
    200-cluster batch, 174 of 414 Refiner packets and 10 MAFFT packets produced
    a consensus of length ZERO, because no column of a 4.2 kb element reached
    50% of ALL rows. Refiner itself had called a consensus in 174/174 of them
    (median 4,217 bp) -- the alignment was fine, the occupancy rule was not.

    A row's span runs from its first to its last non-gap column; interior gaps
    are real deletions and still count against it.
    """
    gap = ord("-")
    nz = m != gap
    n, w = m.shape
    if not n or not w:
        return np.zeros(w), np.zeros(w, dtype=np.int64)
    has = nz.any(axis=1)
    first = np.argmax(nz, axis=1)
    last = w - 1 - np.argmax(nz[:, ::-1], axis=1)
    delta = np.zeros(w + 1, dtype=np.int64)
    np.add.at(delta, first[has], 1)
    np.add.at(delta, last[has] + 1, -1)
    cover = np.cumsum(delta)[:w]
    occ = nz.sum(axis=0) / np.maximum(cover, 1)
    return occ, cover


def consensus(m: np.ndarray, min_occupancy: float = 0.5,
              tie: str = "iupac") -> tuple[str, np.ndarray, np.ndarray]:
    """Majority consensus with an occupancy rule.

    Returns (consensus_string, occupancy_per_column, is_match_column).
    A column with occupancy below min_occupancy is an INSERT column: it is not
    part of the consensus, and in Stockholm terms it is a '.' in the RF line.
    """
    gap = ord("-")
    occ, cover = span_occupancy(m)
    # BOTH tests are needed, and they catch opposite failures. Span-normalised
    # occupancy alone would keep a ragged shoulder that only two rows reach --
    # 2/2 reads as fully occupied. An absolute depth floor alone was what threw
    # away 174 Refiner packets. A match column must be agreed by most rows that
    # reach it AND be reached by enough rows to mean anything.
    depth_floor = min(MIN_MATCH_DEPTH, m.shape[0])
    is_match = (occ >= min_occupancy) & (cover >= depth_floor)
    letters = np.array([ord(c) for c in "ACGT"], dtype=np.uint8)
    counts = np.stack([(m == L).sum(axis=0) for L in letters])   # 4 x width
    best = counts.argmax(axis=0)
    total = counts.sum(axis=0)
    out = []
    for j in range(m.shape[1]):
        if not is_match[j]:
            continue
        if total[j] == 0:
            out.append("N")
            continue
        top = counts[:, j].max()
        if tie == "iupac" and (counts[:, j] == top).sum() > 1:
            out.append(_iupac(counts[:, j] == top))
        else:
            out.append("ACGT"[best[j]])
    return "".join(out), occ, is_match


_IUPAC = {"AG": "R", "CT": "Y", "GT": "K", "AC": "M", "CG": "S", "AT": "W",
          "CGT": "B", "AGT": "D", "ACT": "H", "ACG": "V", "ACGT": "N"}


def _iupac(mask: np.ndarray) -> str:
    key = "".join(b for b, k in zip("ACGT", mask) if k)
    return _IUPAC.get(key, key if len(key) == 1 else "N")


def trim_alignment(m: np.ndarray, min_occupancy: float = 0.5,
                   ) -> tuple[int, int]:
    """Columns [lo, hi) spanned by the well-occupied core.

    Trims the ragged 5'/3' shoulders that a handful of over-long rows create,
    without touching interior low-occupancy columns (those are real deletions).
    Occupancy is span-normalised (see span_occupancy).
    """
    occ, cover = span_occupancy(m)
    depth_floor = min(MIN_MATCH_DEPTH, m.shape[0])
    ok = np.flatnonzero((occ >= min_occupancy) & (cover >= depth_floor))
    if not len(ok):
        return 0, m.shape[1]
    return int(ok[0]), int(ok[-1]) + 1


# ------------------------------------------- 6. keeping ids true after trimming
def _flip(strand: str) -> str:
    return "-" if strand == "+" else "+"


def row_coords_mafft(recs, coords: dict) -> list:
    """Genomic interval + orientation of each MAFFT row AS WRITTEN.

    MAFFT --adjustdirectionaccurately reverse-complements rows it judges
    backwards and marks them '_R_'.  Such a row holds the opposite strand from
    the one its identifier named, so the orientation recorded here flips with
    it (measured: 77 of 100 rows in cluster 540).
    """
    out = []
    for rid, _ in recs:
        flipped = rid.startswith("_R_")
        chrom, s, e, strand = coords[rid[3:] if flipped else rid]
        st = norm_strand(strand)
        out.append((chrom, s, e, _flip(st) if flipped else st))
    return out


def finalize_alignment(recs: list[tuple[str, str]], row_coords: list,
                       min_occupancy: float = 0.5, min_row_bp: int = 30,
                       max_iter: int = 4) -> dict:
    """Trim, drop rows left empty by the trim, and RETAG every identifier.

    `row_coords[i]` is the genomic interval and orientation of row i AS THE ROW
    IS WRITTEN -- callers resolve engine-specific quirks (MAFFT's '_R_' flip,
    Refiner's appended sub-range) before calling.  Here the only adjustment is
    the trim itself: a trimmed row no longer contains the whole interval its
    identifier claims, so the coordinates are walked in by the number of
    non-gap bases actually removed from each end -- from the 3' end first when
    the row is written on the minus strand.

    Without this, `stk lint --genome` rejects the rows; and a seed whose
    identifiers are subtly wrong is worse than one that fails loudly, because
    Dfam's TSD and extension algorithms fetch flanking sequence BY identifier.
    """
    gap = ord("-")
    _, m = msa_matrix(recs)
    keep = np.ones(len(recs), dtype=bool)
    lo, hi = 0, m.shape[1]

    for _ in range(max_iter):
        lo, hi = trim_alignment(m[keep], min_occupancy)
        inside = (m[:, lo:hi] != gap).sum(axis=1)
        new_keep = keep & (inside >= min_row_bp)
        if new_keep.sum() < 2:                 # never trim the alignment away
            new_keep = keep & (inside > 0)
        if (new_keep == keep).all():
            break
        keep = new_keep

    rows_out, coords_out, dropped = [], [], []
    for i, (chrom, s, e, strand) in enumerate(row_coords):
        nl = int((m[i, :lo] != gap).sum())
        nr = int((m[i, hi:] != gap).sum())
        nin = int((m[i, lo:hi] != gap).sum())
        if not keep[i]:
            dropped.append(dict(seq_id=recs[i][0], aligned_bp=nin,
                                reason="empty_after_trim" if nin == 0
                                       else "below_min_row_bp"))
            continue
        if strand == "+":
            coords_out.append((chrom, s + nl, e - nr, "+"))
        else:
            coords_out.append((chrom, s + nr, e - nl, "-"))
        rows_out.append(m[i, lo:hi])

    mat = np.stack(rows_out) if rows_out else np.zeros((0, hi - lo), np.uint8)
    # a column left all-gap by row removal carries nothing, and dropping it
    # moves no bases, so identifiers stay correct
    if len(mat):
        mat = mat[:, (mat != gap).any(axis=0)]
    cons, occ, is_match = consensus(mat, min_occupancy)
    return dict(matrix=mat, coords=coords_out, consensus=cons,
                is_match=is_match, occupancy=occ, dropped=dropped,
                trim_lo=lo, trim_hi=hi)


# ------------------------------------------------------------- 7. Refiner
def refiner(in_fa: Path, refiner_bin: Path, threads: int = 4,
            workdir: Path | None = None) -> dict:
    """Run Dfam's Refiner and return its alignment plus consensus.

    Refiner is NOT standalone despite appearances: it needs RepeatMasker's
    pure-Perl modules on PERL5LIB (not REPEATMASKER_DIR -- `use lib` runs at
    compile time, before RepModelConfig applies environment overrides) and a
    genuine `rmblastn`.  Stock blastn cannot substitute: Refiner's search uses
    -complexity_adjust, -matrix and the kdiv/cpg_kdiv output fields, none of
    which stock BLAST+ supports.  vendor/refiner/ carries a working install and
    tools/setup_refiner.sh reproduces it.

    Refiner writes beside its input, so it is run in a scratch directory.
    Usefully, it PRESERVES the input identifier and appends the sub-range it
    actually used, e.g.
        GCA_...:OX637595.1:3908840-3909106_+:3-267_+
    which is what makes its rows retaggable back to genomic coordinates.
    """
    import shutil
    # resolved: Refiner runs with cwd=work, so a relative binary path would miss
    refiner_bin = Path(refiner_bin).resolve()
    in_fa = Path(in_fa).resolve()
    work = (Path(workdir).resolve() if workdir
            else Path(tempfile.mkdtemp(prefix="refiner_")))
    work.mkdir(parents=True, exist_ok=True)
    local = work / "copies.fa"
    shutil.copy(in_fa, local)
    p = subprocess.run([str(refiner_bin), "-threads", str(threads), str(local)],
                       capture_output=True, text=True, cwd=work)
    stk = local.with_name(local.name + ".refiner.stk")
    cons_f = local.with_name(local.name + ".refiner_cons")
    if p.returncode != 0 or not stk.exists():
        raise RuntimeError(f"Refiner failed ({p.returncode}): "
                           f"{(p.stdout + p.stderr)[-2000:]}")
    ids, rows = [], []
    for line in open(stk):
        line = line.rstrip("\n")
        if not line or line.startswith("#") or line == "//":
            continue
        parts = line.split(None, 1)
        if len(parts) == 2:
            ids.append(parts[0])
            rows.append(parts[1])
    cons = "".join(l.strip() for l in open(cons_f) if not l.startswith(">"))
    kimura = float("nan")
    with open(cons_f) as fh:
        head = fh.readline()
        if "Avg Kimura" in head:
            try:
                kimura = float(head.split("Avg Kimura =")[1].split(")")[0])
            except (IndexError, ValueError):
                pass
    return dict(records=list(zip(ids, rows)), consensus=cons,
                avg_kimura=kimura, workdir=str(work),
                stdout=p.stdout[-500:])


def parse_refiner_id(rid: str, coords: dict) -> tuple | None:
    """'<smitten>:<a>-<b>_<s>' -> genomic (chrom, start, end, strand).

    <a>-<b> is 1-based fully closed ON THE EXTRACTED SEQUENCE, which is the
    genomic forward strand only when the original copy was '+'.  For a '-'
    copy the extracted sequence is already reverse-complemented, so the
    sub-range counts from the genomic 3' end.
    """
    head, _, tail = rid.rpartition(":")
    if not head or "-" not in tail:
        return None
    rng, _, sub_strand = tail.rpartition("_")
    if sub_strand not in ("+", "-") or "-" not in rng:
        return None
    try:
        a, b = (int(x) for x in rng.split("-", 1))
    except ValueError:
        return None
    if head not in coords:
        return None
    chrom, s, e, strand = coords[head]
    if strand == "+":
        ns, ne = s + (a - 1), s + b
    else:
        ns, ne = e - b, e - (a - 1)
    final = strand if sub_strand == "+" else ("-" if strand == "+" else "+")
    return chrom, ns, ne, final


def verify_rows_against_genome(coords: list, rows: list[str],
                               fa: IndexedFasta) -> dict:
    """Check every seed row's sequence against the assembly ourselves.

    `stk lint --genome` reports coordinate problems only as WARN and exits 0,
    and a row whose identifier Smitten cannot parse is dropped from validation
    with NO diagnostic at all -- so a clean lint run is not by itself evidence
    that the coordinates were checked. This does the check directly and counts
    the rows, which makes silence impossible to mistake for success.
    """
    n_ok, bad = 0, []
    for (chrom, s, e, strand), row in zip(coords, rows):
        seq = row.upper().replace(".", "").replace("-", "")
        if chrom not in fa:
            bad.append(dict(chrom=chrom, start=s, end=e, reason="chrom_absent"))
            continue
        ref = fa.fetch(chrom, s, e).upper()
        if strand == "-":
            ref = revcomp(ref)
        if seq == ref:
            n_ok += 1
        else:
            bad.append(dict(chrom=chrom, start=s, end=e, strand=strand,
                            reason="sequence_mismatch",
                            row_len=len(seq), ref_len=len(ref)))
    return dict(n_rows=len(rows), n_verified=n_ok, n_bad=len(bad),
                all_verified=len(bad) == 0, bad=bad[:20])

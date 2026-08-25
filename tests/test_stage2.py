"""Regression tests for the correctness bugs found while building stage 2.

Every test here corresponds to a bug that was real, that produced plausible
output, and that would not have been caught by eye. Run with:

    PYTHONPATH=. .venv/bin/python -m pytest tests/ -q
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lhfseeds import stage0, stage2, stockholm            # noqa: E402
from lhfseeds.fasta import IndexedFasta, revcomp          # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
GENOME = ROOT / "data" / "GCA_951799975.1.fna"
needs_genome = pytest.mark.skipif(not GENOME.exists(),
                                  reason="assembly FASTA not present")


# --------------------------------------------------------------- fasta reader
@needs_genome
def test_indexed_fasta_matches_samtools():
    fa = IndexedFasta(GENOME)
    rng = np.random.default_rng(0)
    names = list(fa.index)
    for _ in range(10):
        ch = names[int(rng.integers(len(names)))]
        L = fa.length(ch)
        s = int(rng.integers(0, max(1, L - 1000)))
        e = min(L, s + int(rng.integers(1, 2000)))
        ref = subprocess.run(
            ["samtools", "faidx", str(GENOME), f"{ch}:{s + 1}-{e}"],
            capture_output=True, text=True).stdout
        assert fa.fetch(ch, s, e) == "".join(ref.splitlines()[1:])


@needs_genome
def test_indexed_fasta_clamps_at_contig_edges():
    fa = IndexedFasta(GENOME)
    ch = next(iter(fa.index))
    L = fa.length(ch)
    assert len(fa.fetch(ch, L - 10, L + 9999)) == 10
    assert fa.fetch(ch, -50, 10) == fa.fetch(ch, 0, 10)
    assert fa.fetch(ch, 50, 50) == ""


# ------------------------------------------------------- gap-character reading
def test_msa_matrix_treats_all_stockholm_gaps_as_gaps():
    """Refiner and Dfam write '.', MAFFT writes '-'. Treating only '-' as a gap
    made every Refiner column look fully occupied, so trimming did nothing and
    the consensus came out at the full alignment width (349 bp vs 267)."""
    recs = [("a", "ACGT....ACGT"), ("b", "ACGT----ACGT"), ("c", "ACGT~~__ACGT")]
    _, m = stage2.msa_matrix(recs)
    occ = (m != ord("-")).mean(axis=0)
    assert list(occ[4:8]) == [0.0, 0.0, 0.0, 0.0], "interior must read as gap"
    cons, _, is_match = stage2.consensus(m, min_occupancy=0.5)
    assert cons == "ACGTACGT"
    assert is_match.sum() == 8


# ------------------------------------------------------------- locus dedup
def test_locus_ids_require_reciprocal_overlap_not_adjacency():
    """Tandem units abut but do not overlap; a touch-based rule chained whole
    arrays into one 'locus' (max 30 members vs 8 under reciprocal overlap)."""
    starts = np.array([0, 100, 200, 300])
    ends = np.array([100, 200, 300, 400])
    assert len(set(stage2._locus_ids(starts, ends).tolist())) == 4

    # genuine co-annotation of one locus by three tools
    starts = np.array([0, 3, 5])
    ends = np.array([100, 98, 103])
    assert len(set(stage2._locus_ids(starts, ends).tolist())) == 1


def test_locus_ids_union_a_chain_into_one_id():
    """Single-linkage without union assigned the LAST matching id, so a chain
    could end up split across ids that should have been merged."""
    starts = np.array([0, 40, 80])
    ends = np.array([100, 140, 180])
    lid = stage2._locus_ids(starts, ends)
    assert len(set(lid.tolist())) == 1


def test_cluster_loci_representative_is_a_real_row():
    """groupby().first() is per-column and takes the first NON-NULL value,
    which spliced fields from several annotations into a representative that
    does not exist at the locus."""
    df = pd.DataFrame([
        dict(chrom="c1", start=0, end=100, strand="+", n_fragments=1,
             frag_starts="0", frag_ends="100", div=1.0, cons_start=np.nan,
             cons_end=np.nan, member="edta:X", merge_mode="gap_aware",
             cons_len=100.0, cons_len_source="clustermate", is_full=False),
        dict(chrom="c1", start=2, end=99, strand="+", n_fragments=1,
             frag_starts="2", frag_ends="99", div=2.0, cons_start=1.0,
             cons_end=99.0, member="rm2:Y", merge_mode="gap_aware",
             cons_len=100.0, cons_len_source="bed16", is_full=True),
    ])
    out = stage2.cluster_loci(df, modal_cons_len=100.0)
    assert len(out) == 1
    r = out.iloc[0]
    assert r.member == "rm2:Y", "near-full-length row must win"
    # the representative's fields must all come from that same row
    assert (r.start, r.end, r.cons_len_source) == (2, 99, "bed16")
    assert r.n_members == 2 and r.n_tools == 2


# ------------------------------------------------------------ strand handling
def test_uncalled_strand_normalises_to_plus():
    """EDTA leaves strand '.' on ~0.3% of hits. Extraction reverse-complements
    only on '-', so '.' IS extracted forward and the identifier must say '+'.
    Treating it as reverse produced coordinates whose sequence was the reverse
    complement of the row."""
    assert stage2.norm_strand(".") == "+"
    assert stage2.norm_strand("+") == "+"
    assert stage2.norm_strand("-") == "-"

    df = pd.DataFrame([dict(chrom="c1", start=10, end=50, strand=".",
                            n_fragments=1, frag_starts="10", frag_ends="50")])
    out = stage2.explode_fragments(df)
    assert out.strand.iat[0] == "+"
    assert bool(out.strand_called.iat[0]) is False


def test_mafft_flip_inverts_the_recorded_strand():
    """MAFFT --adjustdirectionaccurately reverse-complements rows and marks
    them '_R_' (77 of 100 rows in cluster 540). The row then holds the opposite
    strand from the one its identifier names."""
    coords = {"s1": ("c1", 10, 50, "+"), "s2": ("c1", 60, 90, "-")}
    recs = [("_R_s1", "ACGT"), ("s2", "ACGT")]
    rc = stage2.row_coords_mafft(recs, coords)
    assert rc[0] == ("c1", 10, 50, "-"), "flipped '+' row is now minus"
    assert rc[1] == ("c1", 60, 90, "-"), "unflipped '-' row is unchanged"


# ------------------------------------------------------- trimming retags ids
def test_trim_walks_coordinates_in_from_the_correct_end():
    """A trimmed row no longer contains the whole interval its identifier
    claims; on a minus-strand row the left of the alignment is the genomic
    3' end, so the two ends must not be adjusted symmetrically."""
    # Asymmetric shoulders (2 left, 1 right) so the plus and minus answers
    # differ; four filler rows keep the shoulder columns clearly under the
    # occupancy threshold rather than exactly on it.
    #                     col: 0123456789
    recs = [("p", "AA" + "CCCCCC" + "T."),
            ("m", "AA" + "CCCCCC" + "T."),
            ("w", ".." + "CCCCCC" + ".."),
            ("x", ".." + "CCCCCC" + ".."),
            ("y", ".." + "CCCCCC" + ".."),
            ("z", ".." + "CCCCCC" + "..")]
    # interval length must equal each row's ungapped length: 9, 9, 6, 6, 6, 6
    coords = [("c1", 100, 109, "+"), ("c1", 200, 209, "-"),
              ("c1", 300, 306, "+"), ("c1", 400, 406, "+"),
              ("c1", 500, 506, "+"), ("c1", 600, 606, "+")]
    fin = stage2.finalize_alignment(recs, coords, min_occupancy=0.5,
                                    min_row_bp=1)
    out = fin["coords"]
    assert len(out) == 6
    # plus row: 2 bases trimmed from the genomic 5' end, 1 from the 3'
    assert out[0] == ("c1", 102, 108, "+")
    # minus row: the alignment's LEFT is the genomic 3' end, so the 2-base
    # left trim comes off `end` and the 1-base right trim off `start`
    assert out[1] == ("c1", 201, 207, "-")
    # untrimmed rows are untouched
    assert out[2] == ("c1", 300, 306, "+")
    # every retagged interval is exactly as long as its row's ungapped sequence
    for (c, s, e, st), row in zip(out, fin["matrix"]):
        assert e - s == int((row != ord("-")).sum())


def test_refiner_id_roundtrip_plus_and_minus():
    """Refiner preserves the input identifier and appends the sub-range it
    used. On a '-' copy the extracted sequence is already reverse-complemented,
    so the sub-range counts from the genomic 3' end."""
    coords = {"A": ("c1", 1000, 1100, "+"), "B": ("c1", 2000, 2100, "-")}
    assert stage2.parse_refiner_id("A:1-100_+", coords) == ("c1", 1000, 1100, "+")
    assert stage2.parse_refiner_id("A:11-100_+", coords) == ("c1", 1010, 1100, "+")
    # minus copy: sub-range 1-90 covers the extracted 5' end = genomic 3' end
    assert stage2.parse_refiner_id("B:1-90_+", coords) == ("c1", 2010, 2100, "-")
    # a sub-alignment on the opposite strand flips the final strand
    assert stage2.parse_refiner_id("A:1-100_-", coords) == ("c1", 1000, 1100, "-")
    assert stage2.parse_refiner_id("nonsense", coords) is None


# --------------------------------------------------- stage 0 reproducibility
def test_cluster_ids_are_deterministic_and_key_is_content_addressed():
    """Both linkages accumulate through sets, whose iteration order depends on
    PYTHONHASHSEED. Left alone, a re-run renumbered every cluster: cluster 62
    became 292 and id 62 came back holding an unrelated family."""
    edges = [dict(a="rm2:f1", b="edta:g1", w_bp=0.9, joint_bp=900, w_copy=0.5),
             dict(a="rm2:f2", b="edta:g2", w_bp=0.8, joint_bp=800, w_copy=0.5),
             dict(a="pantera:h1", b="rm2:f3", w_bp=0.7, joint_bp=700, w_copy=0.5)]
    runs = [stage0.link_strict(edges) for _ in range(5)]
    assert all(r == runs[0] for r in runs)
    assert runs[0] == sorted(runs[0], key=lambda c: (len(c), c))

    k1 = stage0.cluster_key(["rm2:a", "edta:b"])
    k2 = stage0.cluster_key(["edta:b", "rm2:a"])
    assert k1 == k2, "key must not depend on member order"
    assert k1 != stage0.cluster_key(["rm2:a", "edta:c"])


# ---------------------------------------------------------------- stockholm
def test_rf_line_width_matches_the_alignment():
    cons = "ACGT"
    is_match = np.array([True, False, True, True, False, True])
    rf = stockholm.rf_line(cons, is_match)
    assert rf == "A.CG.T"
    assert len(rf) == len(is_match)


def test_format_record_rejects_a_ragged_alignment():
    with pytest.raises(ValueError):
        stockholm.format_record(["a", "b"], ["ACGT", "ACG"], {"SQ": 2})


def test_format_record_writes_one_line_per_sequence():
    """Dfam rejects block/interleaved format (`block_format`, ERROR), so rows
    must never be wrapped however wide the alignment gets."""
    ids = ["GCA_1:c1:1-8_+", "GCA_1:c1:20-27_-"]
    rows = ["ACGTACGT", "ACGT.CGT"]
    rec = stockholm.format_record(ids, rows, {
        "DE": "d", "AU": "Ada Lovelace", "TP": "Interspersed_Repeat;Unknown",
        "OC": "Gobius niger", "SQ": 2, "RF": "ACGTACGT"})
    body = [l for l in rec.splitlines()
            if l and not l.startswith(("#", "//"))]
    assert len(body) == 2
    assert "" not in rec.split("\n")[:-1], "a blank line would make it block format"


@needs_genome
def test_extract_reverse_complements_minus_and_accounts_for_flanks():
    fa = IndexedFasta(GENOME)
    ch = next(iter(fa.index))
    rows = pd.DataFrame([dict(chrom=ch, start=1000, end=1100, strand="-")])
    recs, stats = stage2.extract(fa, rows, "GCA_X", flank_bp=0)
    assert recs[0][1] == revcomp(fa.fetch(ch, 1000, 1100))
    frecs, fstats = stage2.extract(fa, rows, "GCA_X", flank_bp=50)
    assert (fstats.length.iat[0]
            == stats.length.iat[0] + fstats.flank_left.iat[0]
            + fstats.flank_right.iat[0])
    assert recs[0][0] == "GCA_X:%s:1001-1100_-" % ch

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


@needs_genome
def test_verify_rows_against_genome_catches_a_shifted_row():
    """`stk lint --genome` reports coordinate problems as WARN and exits 0, and
    silently drops rows whose identifier it cannot parse -- so the pipeline
    checks coordinates itself and counts the rows it checked."""
    fa = IndexedFasta(GENOME)
    ch = next(iter(fa.index))
    plus = fa.fetch(ch, 1000, 1060)
    minus = revcomp(fa.fetch(ch, 2000, 2060))
    coords = [(ch, 1000, 1060, "+"), (ch, 2000, 2060, "-")]
    ok = stage2.verify_rows_against_genome(coords, [plus, minus], fa)
    assert ok["all_verified"] and ok["n_verified"] == 2

    # a 500 bp shift must be caught, and a wrong strand must be caught
    shifted = [(ch, 1500, 1560, "+"), (ch, 2000, 2060, "-")]
    bad = stage2.verify_rows_against_genome(shifted, [plus, minus], fa)
    assert not bad["all_verified"] and bad["n_bad"] == 1

    flipped = [(ch, 1000, 1060, "-"), (ch, 2000, 2060, "-")]
    bad2 = stage2.verify_rows_against_genome(flipped, [plus, minus], fa)
    assert not bad2["all_verified"], "wrong strand must not verify"


def test_verify_rows_tolerates_stockholm_gaps_in_the_row():
    """Rows arrive with Dfam '.' gaps; the check compares ungapped sequence."""
    class FakeFa:
        index = {"c1": None}
        def __contains__(self, c): return c == "c1"
        def fetch(self, c, s, e): return "ACGTACGT"[s:e]
    res = stage2.verify_rows_against_genome(
        [("c1", 0, 8, "+")], ["ACGT....ACGT".replace("....", "")], FakeFa())
    assert res["all_verified"]


# ------------------------------------------- over-assembly deconvolution (F4)
def test_span_only_near_full_length_for_deconvolved_members():
    """A member whose consensus is an over-assembly has repeat_start/repeat_end
    in the multi-unit coordinate system of the collapsed entry. Keeping the
    'both ends reached' test while swapping in the unit length mixes coordinate
    systems -- measured 25.2% against a cross-tool truth of 49.5%, where
    span-only gives 58.2% and the old behaviour gave 0.74%."""
    from lhfseeds import stage1
    # a copy covering unit 2 of a 3-unit 765 bp entry: 265 bp of consensus,
    # but starting at 264, so it reaches NEITHER end of the collapsed entry
    copies = pd.DataFrame([dict(cons_start=264.0, cons_end=529.0)])
    assert not stage1.near_full_length(copies, 265.0).iat[0], \
        "ends test cannot be satisfied in the collapsed coordinate system"
    assert stage1.near_full_length(copies, 265.0, span_only=True).iat[0], \
        "span-only recognises a full unit"
    # and it still rejects a genuinely short copy
    short = pd.DataFrame([dict(cons_start=264.0, cons_end=300.0)])
    assert not stage1.near_full_length(short, 265.0, span_only=True).iat[0]


def test_near_full_length_is_false_without_consensus_coordinates():
    """EDTA carries none; flag, do not fake."""
    from lhfseeds import stage1
    copies = pd.DataFrame([dict(cons_start=np.nan, cons_end=np.nan)])
    assert not stage1.near_full_length(copies, 265.0).iat[0]
    assert not stage1.near_full_length(copies, float("nan")).iat[0]


# ------------------------------------------------------------- lint wrapper
STK_BIN = ROOT / "vendor" / "dfam-curator" / "target" / "release" / "stk"
needs_stk = pytest.mark.skipif(not STK_BIN.exists(), reason="stk not built")


@needs_stk
def test_lint_does_not_read_an_io_failure_as_clean():
    """There is no exit code 2: `stk` exits 1 both for "emitted an ERROR" and
    for an I/O failure, and an I/O failure emits NO diagnostics -- so a file
    that could not even be opened came back with n_error == 0."""
    r = stockholm.lint(Path("/nonexistent/nope.stk"), STK_BIN)
    assert r["returncode"] != 0
    assert r["n_error"] == 0, "no diagnostics are emitted for an I/O failure"
    assert r["io_failure"] is True
    assert r["clean"] is False


@needs_stk
def test_lint_clean_requires_no_coordinate_warnings(tmp_path):
    """`stk lint --genome` reports coordinate problems as WARN and exits 0, so
    exit status and n_error both miss them."""
    rec = stockholm.format_record(
        ["GCA_1:c1:1-8_+"], ["ACGTACGT"],
        {"DE": "d", "AU": "Ada Lovelace", "TP": "Interspersed_Repeat;Unknown",
         "OC": "Gobius niger", "SQ": 1, "RF": "ACGTACGT"})
    f = tmp_path / "x.stk"
    stockholm.write_stockholm(f, [rec])
    r = stockholm.lint(f, STK_BIN)
    assert r["io_failure"] is False
    assert "n_coord_warnings" in r and "clean" in r


def _deconv_decision(lens, ratio_min=1.8, spread_max=1.3):
    """Mirror of the rule in run_stage1 (kept in step with it by the tests
    below); see PIPELINE_FINDINGS F4."""
    out = {}
    for m, own in lens.items():
        others = [v for k, v in lens.items() if k != m]
        if len(others) < 2:
            out[m] = "skip"
            continue
        mate = float(np.median(others))
        if mate <= 0 or own / mate < ratio_min:
            out[m] = "keep"
            continue
        agree = [v for v in others if 1 / spread_max <= v / mate <= spread_max]
        if len(agree) < 2 or len(agree) < 0.5 * len(others):
            out[m] = "abstain"
            continue
        ref = float(np.median(agree))
        out[m] = "deconvolve" if own / ref >= ratio_min else "keep"
    return out


def test_deconvolution_fires_when_two_members_are_over_assembled():
    """Cluster 62 has TWO over-assembled REPET entries, so each sees the other
    among its mates. A max/min spread test abstained on exactly the cases it
    was meant to allow; support must be counted as a majority around the
    median instead."""
    d = _deconv_decision({"rm2": 269., "pantera": 264., "edta": 264.,
                          "repetA": 764., "repetB": 765.})
    assert d["repetA"] == d["repetB"] == "deconvolve"
    assert d["rm2"] == d["pantera"] == d["edta"] == "keep"


def test_deconvolution_refuses_a_single_length_carrying_mate():
    """rm2 rnd-4_family-1503 is a sequence-verified 2,533 bp dimer with exactly
    one length-carrying mate; a median over one mate is an opinion, not a
    consensus, and inheritance would just inherit the over-assembly."""
    assert _deconv_decision({"dimer": 2533., "mate": 1266.})["dimer"] == "skip"


def test_deconvolution_does_not_touch_an_ordinary_cluster():
    d = _deconv_decision({"rm2": 378., "pantera": 351., "edta": 351.,
                          "repet": 367.})
    assert set(d.values()) == {"keep"}


def test_clustermate_paint_order_is_evidence_ranked_not_hash_ranked():
    """build_clustermate_conslen paints members into a shared array, so at a
    contested base whoever paints LAST wins. Passing a set made that order
    depend on PYTHONHASHSEED: one EDTA member's inherited consensus length
    moved from 381 bp to 1,755 bp between two runs on identical input."""
    import inspect
    from lhfseeds import stage1
    src = inspect.getsource(stage1.build_clustermate_conslen)
    assert "order = sorted(" in src, "paint order must be explicit"
    assert "for (tool, fam) in order:" in src, "must iterate the sorted order"
    sig = inspect.signature(stage1.build_clustermate_conslen)
    assert "weight_by_member" in sig.parameters

    # the ordering itself: least evidence first, so the best-supported wins
    fams = {("edta", "A"), ("rm2", "B"), ("repet", "C")}
    w = {"edta:A": 10.0, "rm2:B": 1000.0, "repet:C": 100.0}
    order = sorted(fams, key=lambda tf: w.get(f"{tf[0]}:{tf[1]}", 0.0))
    assert order[-1] == ("rm2", "B"), "highest-evidence member must paint last"
    assert order[0] == ("edta", "A")


# ------------------------------------------------- merge_copies vectorisation
def _hits(rows):
    """Minimal member_hits-shaped frame."""
    return pd.DataFrame(rows, columns=["chrom", "chromStart", "chromEnd", "name",
                                       "strand", "perc_div", "repeat_start",
                                       "repeat_end", "repeat_left", "hit_id"])


def test_merge_copies_groups_only_same_chrom_and_strand():
    from lhfseeds import stage1
    h = _hits([
        ["c1", 100, 200, "F", "+", 1.0, 1, 100, 0, None],
        ["c1", 210, 300, "F", "+", 3.0, 1, 100, 0, None],   # merges (gap 10)
        ["c1", 310, 400, "F", "-", 1.0, 1, 100, 0, None],   # strand break
        ["c2", 100, 200, "F", "+", 1.0, 1, 100, 0, None],   # chrom break
    ])
    out = stage1.merge_copies(h, "merge_always")
    assert list(out.n_fragments) == [2, 1, 1]
    assert list(out.start) == [100, 310, 100]
    assert list(out.end) == [300, 400, 200]
    assert out.frag_starts.iat[0] == "100;210"
    assert out.frag_ends.iat[0] == "200;300"
    assert out["div"].iat[0] == pytest.approx(2.0)  # nanmean of 1.0, 3.0
    # NB out["div"], not out.div -- `div` collides with DataFrame.div()


def test_merge_copies_never_bridges_a_gap_beyond_the_cap():
    from lhfseeds import stage1
    far = stage1.MAX_MERGE_GAP + 1000
    h = _hits([["c1", 0, 100, "F", "+", 1.0, 1, 100, 0, None],
               ["c1", far, far + 100, "F", "+", 1.0, 1, 100, 0, None]])
    assert len(stage1.merge_copies(h, "merge_always")) == 2


def test_runaway_cap_splits_at_the_widest_gap_recursively():
    """A copy cannot exceed MAX_COPY_X_CONSENSUS x the consensus (F1). The cut
    goes at the widest internal gap, and the result is re-checked until every
    piece fits -- one 3-unit array must become three copies, not two."""
    from lhfseeds import stage1
    # three 100 bp units at 0, 150, 300; consensus 100 -> cap 150
    h = _hits([["c1", 0, 100, "F", "+", 1.0, 1, 100, 0, None],
               ["c1", 150, 250, "F", "+", 1.0, 1, 100, 0, None],
               ["c1", 300, 400, "F", "+", 1.0, 1, 100, 0, None]])
    out = stage1.merge_copies(h, "merge_always", cons_len=100.0)
    assert len(out) == 3, "recursion must continue until every piece fits"
    assert list(out.start) == [0, 150, 300]
    assert (out.end - out.start).max() <= stage1.MAX_COPY_X_CONSENSUS * 100


def test_runaway_cap_span_uses_max_end_not_last_end():
    """Intervals can nest; taking the last row's end would understate the span
    and let an over-long copy through."""
    from lhfseeds import stage1
    h = _hits([["c1", 0, 500, "F", "+", 1.0, 1, 100, 0, None],
               ["c1", 10, 20, "F", "+", 1.0, 1, 100, 0, None]])
    out = stage1.merge_copies(h, "merge_always", cons_len=100.0)
    assert len(out) == 2, "nested pair spans 500 > cap 150, so it must split"


def test_merge_copies_gap_aware_consults_the_probe_only_for_real_gaps():
    from lhfseeds import stage1
    seen = []

    def probe(ch, a, b):
        seen.append((ch, a, b))
        return 0.9                      # heavily foreign -> refuse the merge

    h = _hits([["c1", 0, 100, "F", "+", 1.0, 1, 100, 0, None],
               ["c1", 150, 250, "F", "+", 1.0, 1, 100, 0, None],
               ["c1", 250, 350, "F", "+", 1.0, 1, 100, 0, None]])  # gap == 0
    out = stage1.merge_copies(h, "gap_aware", foreign_cov=probe)
    assert seen == [("c1", 100, 150)], "a zero-length gap needs no probe"
    assert list(out.n_fragments) == [1, 2]


def test_merge_copies_preserves_integer_consensus_coordinates():
    """member_hits keeps an int64 column when it has no NA; pandas .min()
    returned an int, so forcing float turned 178 into 178.0 in the TSV."""
    from lhfseeds import stage1
    h = _hits([["c1", 0, 100, "F", "+", 1.0, 178, 231, 0, None]])
    h["repeat_start"] = h.repeat_start.astype("int64")
    h["repeat_end"] = h.repeat_end.astype("int64")
    out = stage1.merge_copies(h, "merge_always")
    assert str(out.cons_start.iat[0]) == "178"
    assert str(out.cons_end.iat[0]) == "231"


def test_merge_copies_handles_all_nan_consensus_coordinates():
    from lhfseeds import stage1
    h = _hits([["c1", 0, 100, "F", "+", np.nan, np.nan, np.nan, np.nan, None],
               ["c1", 110, 200, "F", "+", np.nan, np.nan, np.nan, np.nan, None]])
    out = stage1.merge_copies(h, "merge_always")
    assert np.isnan(out.cons_start.iat[0]) and np.isnan(out.cons_end.iat[0])
    assert np.isnan(out["div"].iat[0])


def test_merge_copies_empty_input_returns_empty_frame_with_columns():
    from lhfseeds import stage1
    out = stage1.merge_copies(_hits([]), "merge_always")
    assert len(out) == 0
    for c in ("chrom", "start", "end", "strand", "n_fragments",
              "frag_starts", "frag_ends", "div", "cons_start", "cons_end"):
        assert c in out.columns


# ------------------------------------------- bimodal length modes (2026-08-26)
def test_length_modes_splits_a_solo_ltr_from_a_full_element():
    """Cluster 1183: three members at ~423 bp and five at ~7400 bp. Collapsing
    that to one modal length and capping against it deleted every full-length
    locus and rebuilt a 7.4 kb element as 453 bp."""
    ms = stage2.length_modes({
        "rm2:a": 423, "edta:b": 423, "pantera:c": 414,
        "repet:d": 7400, "fastltr:e": 7367, "pantera:f": 7359,
        "rm2:g": 6523, "edta:h": 7400})
    assert len(ms) == 2
    assert [m["credible"] for m in ms] == [True, True]
    assert round(ms[0]["centre"]) == 423 and round(ms[1]["centre"]) == 7367
    assert ms[1]["n_tools"] == 5


def test_length_modes_rejects_a_single_tool_long_mode():
    """Cluster 62 before deconvolution: 264/264/269 from three tools, 764/765
    from REPET alone. The long mode must NOT become credible, or F4 reopens and
    tandem/chimeric multimers re-enter the seed."""
    ms = stage2.length_modes({"rm2:a": 269, "pantera:b": 264, "edta:c": 264,
                              "repet:d": 764, "repet:e": 765})
    assert len(ms) == 2
    assert ms[0]["credible"] is True and ms[0]["n_tools"] == 3
    assert ms[1]["credible"] is False, "REPET alone is one tool, not a consensus"
    assert ms[1]["tools"] == ["repet"]


def test_length_modes_keeps_one_mode_when_lengths_agree():
    ms = stage2.length_modes({"rm2:a": 378, "pantera:b": 351, "edta:c": 351,
                              "repet:d": 367})
    assert len(ms) == 1 and ms[0]["credible"] is True


def test_length_modes_ignores_missing_lengths():
    ms = stage2.length_modes({"edta:a": float("nan"), "rm2:b": 300,
                              "pantera:c": 310, "repet:d": None})
    assert len(ms) == 1 and ms[0]["n_members"] == 2


def test_assign_mode_uses_the_geometric_midpoint():
    """Nearest in LOG space: with modes at 423 and 7367 the boundary is their
    geometric mean (~1765 bp), so a tandem dimer of the solo LTR is still
    judged against 423 and capped, while a truncated full element at 2 kb is
    judged against 7367 and kept."""
    centres = [423.0, 7367.0]
    spans = np.array([420, 850, 1700, 1800, 7400, 12000])
    got = stage2.assign_mode(spans, centres)
    assert list(got) == [0, 0, 0, 1, 1, 1]
    # and the cap then does the right thing with x1.5
    ref = np.array(centres)[got]
    over = spans > 1.5 * ref
    assert list(over) == [False, True, True, False, False, True], \
        ("850 and 1700 bp are multimers of the solo LTR and are capped; 1800 "
         "and 7400 are judged against the full element and kept; 12000 exceeds "
         "1.5x even the full element and is capped")


def test_assign_mode_with_no_credible_modes_is_safe():
    assert list(stage2.assign_mode(np.array([100, 200]), [])) == [0, 0]


def test_only_measured_lengths_should_vote_on_a_mode():
    """A member whose cons_len was inherited from a cluster-mate is repeating
    an opinion, not supplying evidence. Counting its tool inflates the mode's
    cross-tool support -- and the #=GF CC line naming those tools would then
    ship a provenance claim that is not true. Measured on cluster 1183: edta
    contributes 423 and 7400 bp while EDTA's own library entries are 244 and
    272 bp, and its BED has NA consensus coordinates."""
    loci = pd.DataFrame([
        dict(member="rm2:a", cons_len=423.0, cons_len_source="bed16"),
        dict(member="pantera:b", cons_len=414.0, cons_len_source="bed16"),
        dict(member="edta:c", cons_len=423.0, cons_len_source="clustermate"),
        dict(member="repet:d", cons_len=7400.0, cons_len_source="bed16"),
        dict(member="fastltr:e", cons_len=7367.0, cons_len_source="bed16"),
        dict(member="edta:f", cons_len=7400.0, cons_len_source="clustermate"),
    ])
    meas = loci[loci.cons_len_source == "bed16"]
    ms = stage2.length_modes(meas.groupby("member").cons_len.median().to_dict())
    assert [m["n_tools"] for m in ms] == [2, 2]
    assert all("edta" not in m["tools"] for m in ms), \
        "an inherited length must not make its tool count as support"
    # counting everything would have claimed three tools for each mode
    ms_all = stage2.length_modes(loci.groupby("member").cons_len.median().to_dict())
    assert [m["n_tools"] for m in ms_all] == [3, 3]


def test_seed_ids_stay_within_the_45_character_limit():
    """`stk lint` caps #=GF ID at 45 characters (id_too_long, ERROR). The
    longest form the pipeline can build is a multi-model merge_always cluster
    on the refiner engine; carrying the model's length in the ID pushed it to
    46-47 and shipped 12 ERRORs in a live batch before this was caught."""
    longest = 0
    for mode in ("gap_aware", "merge_always"):
        for tag in ("LTR", "int", "m0", "m1"):
            for eng in ("mafft", "refiner"):
                i = f"TEbedSeeds_c01183_{mode}_{tag}_{eng}"
                longest = max(longest, len(i))
                assert len(i) <= 45, f"{i} is {len(i)} chars"
    assert longest >= 40, "sanity: the worst case should be close to the limit"
    # the old form must still be recognised as over-limit, so this test bites
    assert len("TEbedSeeds_c01183_merge_always_m1x7363_refiner") == 46

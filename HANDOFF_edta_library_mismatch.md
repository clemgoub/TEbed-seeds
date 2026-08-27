# Handoff — EDTA family names do not join its library to its annotation (2026-08-26)

> ## UPDATE (same day) — the track is CURRENT; the mismatch is inside EDTA's own output
>
> CG pointed at the GenomeArk source
> (`downstream_analyses/repeats/systematic_annotations/EDTA-v2.3.2/GCA_951799975.1/`)
> and asked: if the library there does not fit the track, the track is stale.
> Checked. **The track is not stale — it reproduces the published annotation**
> for 4 of 5 tested families exactly (TE_00002414 724/724, TE_00001838 30/30,
> TE_00000020 62/62, TE_00000077 80/80, identical longest lengths). One family
> has a small residual: **TE_00000236 has 665 published loci vs 660 in our
> track, longest 4,414 vs 4,404 bp** — cause not established (ingest filtering,
> or the v2.3.1/v2.3.2 discrepancy noted in point 6). Worth a look, but far too
> small to explain the library mismatch. The conclusion below about our local
> `libs/` copy stands, and the root cause is upstream and worse:
>
> 1. The GenomeArk `TElib.fa` (3,389 entries) is a **different file** from our
>    `libs/` copy (6,862 entries): 2,964 shared names, **0 identical
>    sequences**. Our copy should be replaced with the GenomeArk one regardless.
> 2. The GenomeArk library is genuinely for this assembly — all 8 tested
>    entries blast to it genome-wide, 7 of them at **98.9–100%**; the eighth
>    (`TE_00001583`) tops out at **80.9% over 492 of its 570 bp**, so the range
>    is 80.9–100%.
> 3. But it still only matches **15/40** of its own annotated loci (vs 0/40 for
>    our copy). Splitting by EDTA's `method` field does not rescue it:
>    homology loci 2/15, structural loci 0/8.
> 4. **Why:** resolving library entries by coordinate shows the names are
>    numbered independently in the two files. Library `TE_00001838` is
>    annotation `TE_00001858`; library `TE_00002414` is `TE_00002426`. Observed
>    offsets: +0, +8, +10, +11, +12, +20, **+1162** — not a constant shift.
>    Name sets differ too (466 annotated families have no library entry).
> 5. Related, and visible without any library: **47.7% of family names in the
>    published GFF3 carry more than one `classification`** — typically a few
>    `method=structural` loci with one class and many `method=homology` loci
>    with another (e.g. `TE_00000236`: 659 homology `DNA/Helitron` + 6
>    structural `LTR/Gypsy`). This is the real source of the `MITE/DTA`-on-LTR
>    puzzle in F15, and it means **`name` is not a family key in EDTA output**.
> 6. Note the GFF3 header says `EDTA v2.3.1` while the S3 path says
>    `EDTA-v2.3.2`, and its date is 2026-06-24 (objects uploaded 2026-07-31).
>
> **Action:** replace `libs/GCA_951799975.1.fa.mod.EDTA.TElib.fa` with the
> GenomeArk copy (done: `genomeark_edta/` in the deck session), and treat EDTA
> consensus sequences as **unjoinable to EDTA annotations by name**. To attach a
> sequence to an EDTA locus, go through coordinates (extract from the assembly)
> rather than the library. Worth raising with the EDTA authors — this is a
> reproducibility issue in the published output, not a VGP packaging error.
> Evidence: `edta_genomeark_library_test.tsv`, `edta_name_sets.tsv`.

# Original note — the EDTA library on disk does not correspond to edta.bed (2026-08-26)

**Prompted by CG:** *"why EDTA models don't show in the TE-Aid plot, I thought
you had the library"*. We do have a file. It is the wrong file — or rather, it
is a library from a different EDTA run than the one that produced the
annotations. Same failure class as the pantera case in the F14 retraction, and
this time it is measurable at scale.

## Measured

`libs/GCA_951799975.1.fa.mod.EDTA.TElib.fa` (6,862 entries; 2,934 of the 3,791
family names in `inputs/edta.bed` are present by name).

| test | result |
|---|---|
| EDTA library entry vs GENOMIC sequence at that family's own best-length-matched locus | **0 / 40** random families produce any blastn hit (evalue 1e-5, dust off) |
| rm2 library vs its own loci (control) | 9 / 9, pident 91.8–98.1% |
| REPET library vs its own loci (control) | 6 / 6, pident 89.2–96.2% |
| EDTA **LTRlib** coordinate-named entries vs their own loci | 5 / 5 at **100%** — those entries literally are genomic sequence |

The names match; the sequences do not. Several cases match in length exactly
(285 vs 285 bp, 731 vs 731 bp) yet share no homology — the signature of a name
collision between two runs, not of a biological mismatch.

Secondary symptom, visible without the assembly: **52.6% of EDTA families have
a locus longer than 3x their own library consensus** (23% exceed 10x), vs 5.0%
for rm2 and ~0% for REPET and fastltr. It concentrates in the short entries
(74.5% of `MITE/DTA` families), and for 416 of 681 checked families the BED
class disagrees with the library header (e.g. `TE_00000284` annotated
`LTR/Gypsy`, library header `DNA/DTA`). Table:
`edta_lib_vs_annotation_mismatch.tsv`.

## Consequences

1. **Every EDTA-derived sequence comparison in this project is void**, not
   negative. In the rebuilt-vs-member identity test the EDTA members score
   0/20 while rm2 scores 24/26 and REPET 24/25 — that zero is the library, not
   the rebuild. Do not report it as EDTA disagreeing with the consensus.
2. TE-Aid cannot render EDTA member sheets (a 244 bp consensus cannot display
   a 7.4 kb locus), which is why they are absent from the curator demo.
3. **F15 needs revisiting.** F15 reads the `MITE/DTA` header on cluster 1183's
   7.4 kb LTR loci as EDTA reusing a family name. That may be right, but it is
   not established by these data: if the library is from another run, the
   header belongs to a different sequence entirely and says nothing about the
   annotations. Re-derive F15 after the correct library is in hand.
4. F4/C3 per-tool sequence rates: same caveat as pantera — the library that
   produced the BED is required.

## What is needed

The `*.TElib.fa` from **the same EDTA run** that produced
`inputs/edta.bed` (and the same for pantera). Until then, treat EDTA and
pantera as coordinate-only participants: they can join clusters, vote on
support, and contribute copies, but no claim about their consensus sequences is
supportable.

Worth checking whether the hub's EDTA ingest recorded the source run
(`VGP_TEbed/docs/` per-tool provenance) — if the BED came from GenomeArk's
`systematic_annotations/`, the matching library should be retrievable from the
same path.

Data: `edta_library_provenance_evidence.tsv`,
`edta_lib_vs_annotation_mismatch.tsv`, `rebuilt_vs_member_identity.tsv`.

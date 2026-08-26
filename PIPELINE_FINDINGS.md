# Pipeline findings — measured, destined for the slides

Running record of findings from building the LHF pipeline that belong in
`LHF_stage_algorithms.pptx` and (later) the VGP report's technical level.
**Each entry: what was measured, on what, and what it implies.** Update the
deck from this file; re-measure the numbers when the pipeline changes.

Assembly throughout: GCA_951799975.1 (black goby). Tools: the 5
clustering-capable programs (rm2, edta, pantera, fastltr, repet).

---

## F1. Runaway merge chains — transitive merging swallows tandem arrays

**Measured (stage 1, cluster 62):** merging adjacent same-family hits across
gaps produced single "copies" up to **329× the family's consensus length**
(88,525 bp against a 269 bp consensus), with **32% of merge-always copies
longer than 1.5× the consensus**.

**Fix:** cap a copy at `MAX_COPY_X_CONSENSUS` (1.5) × consensus length,
splitting at the widest internal gap until every piece is within bound.
Max copy 88.5 kb → **14.4 kb**; fraction over 1.5× → 0 with a real cap in
force.

**Why it matters for the slides:** this is the copy-level face of a general
hazard — *transitive grouping has no natural stopping point in a repeat-dense
genome*. It is the same failure mode as single-linkage chaining in the
clustering graph (F2), one level down. Both need an explicit bound.

**Slide home:** stage 1.1 (`d5_fill.png` / merging), as a second panel or a
notes item beside the fill measurement.

---

## F2. Chaining in the clustering graph (already in the deck)

Connected components with REPET included: largest cluster **125 members**,
**33%** of clusters order-mixed. Strict (mutual-best) linkage: max **14**,
density 1.00, 25% order-mixed. Lenient (density-split): max 76, density 0.81.

Caveat recorded per user: order-mixing assumes tool classifiers are correct,
which is not fully true. It is used **comparatively** (all strategies scored
against the same labels, so classifier error cancels in the ranking) and is
advisory per cluster, never a filter.

**Slide home:** stage 0.3 (`d3_linkage.png`) — already there.

---

## F3. Gap-aware merging costs almost no usable evidence

**Measured (clusters 62 and 540):**

| cluster | mode | copies | near-full-length | multi-fragment copies |
|---|---|---|---|---|
| 62 | merge_always | 83,692 | 29,136 | 881 |
| 62 | gap_aware | 84,369 | 28,845 | **214** |
| 540 | merge_always | 22,762 | 6,870 | 1,060 |
| 540 | gap_aware | 23,633 | 6,515 | **202** |

Gap-aware yields **more, shorter copies** and cuts multi-fragment copies
**4–5×**, while near-full-length counts barely move (29,136 vs 28,845, −1%).

**Implication:** declining to bridge foreign-filled gaps costs ~1% of the
full-length evidence and removes most of the suspicious multi-fragment
constructions. Provisional default: **gap_aware**.

**Confirmed at seed level — see F8.** The comparison expected to be decisive
(seed depth and rebuilt consensus length) turned out **not** to separate the
modes at all: identical consensus length and depth. Gap-aware is kept as the
default on the secondary criteria instead.

**Slide home:** stage 1.1, paired with the 62.9% foreign-fill measurement.

---

## F4. A tool consensus can be an over-assembly — CHIMERA, not tandem (corrected 2026-08-25); confirmed for four REPET entries, and NOT generalisable

**Prompted by:** cluster 62's two REPET families having suspiciously clean
consensus lengths of **764 and 765 bp**, against 264–269 bp from rm2, pantera
and edta for the same element, and the hypothesis that this was a
palindromic/MITE structure worth tolerating.

### CORRECTION (2026-08-25): the mechanism is a CHIMERA of two elements

The tandem reading below was wrong, and was caught by looking at the locus in
the genome browser. The earlier analysis filtered each tool to the single
cluster family, which hid everything those tools call under *other* family
names at the same locus.

**The decisive test** (co-annotation identity per consensus block; tandem units
of one element must be co-annotated by the SAME family at every block, a
chimera by different families):

| REPET consensus block | dominant rm2 co-annotation | purity |
|---|---|---|
| `G1473-Map8` 51-484 | `rnd-1_family-608` (**496 bp** consensus) | 55-100% |
| `G1473-Map8` 484-765 | `rnd-1_family-286` (**269 bp** consensus) | 79-100% |
| `G1303-Map20` 1-254 | `rnd-1_family-286` (269 bp) | 99-100% |
| `G1303-Map20` 280-764 | `rnd-1_family-608` (496 bp) | 51-100% |

**496 + 269 = 765.** The junction falls exactly where the arithmetic predicts,
and the two REPET entries are the *same pair fused in opposite orders* — which
explains why they partition loci almost perfectly rather than double-annotating.

Both parts are genuine, independently confirmed elements, each co-annotated
along its whole length by a matching-length family from another tool:

- **element A ≈ 496 bp** — rm2 `rnd-1_family-608` (median copy 473 bp = 95% of
  its consensus) and pantera `Unknown_390-fGobNig` (497 bp consensus, median
  copy 482 bp = 97%), which covers A end to end.
- **element B ≈ 269 bp** — rm2 `rnd-1_family-286` (median copy 224 bp = 83%)
  and pantera `Unknown_572-fGobNig` (264 bp, median 237 bp = 90%).

Worked example, matching the browser: `OX637613.1:25,978,734-25,979,682`
(UCSC `chr19:25,978,734-25,979,682`). One REPET copy (`hit_id ms119910`,
consensus 2-764 of 765) spans **both** elements: its 497 bp fragment sits on A
(where pantera calls `Unknown_390` at 100% of its 497 bp consensus and edta
calls `TE_00001450`) and its 277 bp fragment sits on B (rm2 `rnd-1_family-286`
at 99% of 269 bp, pantera `Unknown_572` at 100% of 264 bp, edta `TE_00002952`).
Table and figure: `runs/<asm>/examples/overassembly_locus_example.tsv`.

**Why the coordinate-only phase signal was fooled:** 496 ≈ 2 x 250, so
endpoints cluster near multiples of ~255 and a k=3 array is mimicked by a
chimera whose two parts happen to be ~2x and ~1x a similar length. Phase
enrichment cannot separate those; co-annotation identity can.

**Consequence for the deconvolution rule (§4.1).** Substituting the mate length
is still the right action, but it is right for only *part* of a chimeric
family's copies. Cluster 1103 contains the B-side mates (264, 269) and not the
A-side ones, so REPET copies get judged against ~269 — correct for the B
fragments, wrong for the ~21% of `G1473-Map8` hits (683 of 3,280) that fall in
the A-dominated blocks. The rule should therefore be reported as *coarse* and
the per-copy `cons_len_source` retained, which it is.

### SUPERSEDED reading (kept for the record): ~3 tandem units of a ~255 bp element

- **Palindromic structure: refuted.** Adjacent rm2 hits within 700 bp are
  **75% same-strand**; a palindrome would show an *excess* of opposite-strand
  pairs, not a deficit.  (This part stands — a palindrome is ruled out either way.)
- **rm2's 269 bp entry is the UNIT, not a fragment.** Self-alignment shows no
  same-strand off-diagonal at any offset (tandem cover 0.00 for k=2…5). Its
  only off-diagonal is a **14 bp opposite-strand terminal pair — a terminal
  inverted repeat** (10 bp perfect, `CTTTAAAGGG`). It is a TIR element / MITE.
- **k=3 is uniquely identified for the REPET entries.** Phase enrichment for
  `G1303-Map20` (L=764) across k = 2 / 2.5 / 3 / 3.5 / 4 / 4.5 / 5 / 6 is
  0.053 / 1.152 / **2.839** / 0.245 / 0.569 / 0.449 / 0.737 / 0.854 — k=3 is
  the only value above 1.16. `G1473-Map8` (L=765) behaves the same. Their
  internal endpoint modes sit at 260 bp and 500 bp: the two boundaries of a
  3 × 255 bp array.
- **Two further instances found**, not previously known:
  `Gnig_TEdenovoGr-B-G1039-Map3` (765 bp) and `Gnig_TEdenovoGr-B-G1568-Map3`
  (771 bp), in two other clusters, both k=3, unit 254.7–257.0 bp, each with
  rm2/pantera mates at 264–270 bp. Four REPET entries in total.

### CORRECTION (2026-08-25, REPET library in hand): it is not a tandem array — it is a CHIMERA of two different elements

The sequence test §C was waiting for now runs, and it **overturns the mechanism**
while leaving the pipeline consequence intact. Three independent lines:

1. **Self-alignment finds no tandem structure.** `G1303-Map20` (764 bp) against
   itself (`blastn -word_size 7 -dust no -strand both`) has **no same-strand
   off-diagonal at ±255 or ±510** — the predicted signature of three tandem
   units. Its longest off-diagonal is 22 bp. The only notable hits are
   *opposite*-strand terminal pairs (`2-13` vs `267-256`, 12 bp, 91.7%) — the
   TIR of ONE ~267 bp unit.
2. **The unit matches the long entry exactly ONCE.** rm2's 269 bp
   `rnd-1_family-286` aligns to `G1303-Map20` at **positions 2–269, 98.5% over
   268 bp — and nowhere else**. Against the sibling entry `G1473-Map8` it aligns
   once at **positions 497–764, 96.6%**. One copy each, at opposite ends.
3. **The other ~496 bp is a different, real TE — a second rm2 family.**
   `G1303-Map20`'s positions 270–764 match rm2 **`rnd-1_family-608`** at
   **98.8% over 492 bp** (best hit by bitscore), and a near-identical relative
   `rnd-5_family-6390#SINE/5S-Deu-L2` at 97.4%, so the cassette is a
   5S-derived SINE. The arithmetic is exact: **496 + 269 = 765**. The same
   cassette appears in at least six further REPET entries (`G1329-Map5`,
   `G339-Map20`, `G1058-Map5`, `G926-Map14`, `G652-Map5`) at 99%+ identity and
   at varying offsets. And `G1473-Map8`'s non-element half is the **reverse
   complement** of `G1303-Map20`'s (98.99% over 496 bp).

   *Independently confirmed by CG in the genome browser* at
   `OX637613.1:25,978,734-25,979,682`: the two REPET entries are
   `rnd-1_family-608` + `rnd-1_family-286` fused, in opposite orders.

**What REPET actually did:** a ~265 bp TIR MITE and a ~495 bp 5S-derived SINE
insert next to each other frequently, and REPET built library entries that
**fuse the two**, in both relative orientations — `G1303` = MITE+SINE,
`G1473` = SINE(rc)+MITE. That is a *chimeric consensus*, not a collapsed tandem
array, and it explains every coordinate observation that motivated the tandem
reading: copies cover only the MITE third, the endpoint modes at ~260 and ~500
are the fusion boundary, `L/median_copy ≈ 3` is arithmetic rather than a repeat
count, and the two entries partition loci because they are two constructions of
the same pair.

**So `k` was never real.** The phase-enrichment peak at k=3 was detecting the
fusion boundary, not a period. This is consistent with — and explains — the
decoy-null result (no excess of integer-k structure, p = 0.72) and the
0/8 sensitivity on sequence-verified rm2 cases: the detector was built for a
phenomenon that is not the one occurring. It also retrospectively justifies the
decision **not** to divide by a detected k.

**The pipeline consequence is unchanged and if anything better founded:** the
member's own consensus length must not judge near-full-length, because only
~265 of its 764 bp is the element. Inheriting the cluster-mates' length is
right for a chimera exactly as it would have been for a tandem array. The
provenance value `cons_len_source = deconvolved` keeps its meaning ("this
member's length was overridden by cross-tool consensus"); only the mechanism
named in the docs changes.

**Chimeras and tandem arrays must be counted separately, and only sequence can
separate them.** The coordinate-only phase detector **conflates** the two: a
496 bp cassette fused to a 269 bp element gives L/median_copy ≈ 3 and an
endpoint mode near L/3, which is indistinguishable from a genuine 3-unit array
on coordinates alone. That, not merely low power, is why the decoy null came
back at p = 0.72.

**The real §C3 detector:** a library consensus whose segments align end-to-end
to **two or more different shorter consensi covering disjoint parts of its
length** is a chimera; internal same-strand repeat structure is a tandem
over-assembly. That is a sequence test over the libraries,
it needs no coordinates, and unlike the coordinate screen it targets the
phenomenon that actually occurs. Four of the five clusterable tools now have a
goby library (rm2, edta, repet, fastltr — pantera is still missing), so a
per-tool chimera rate is finally computable.

### Refuted: there is no measurable per-tool over-assembly rate

A coordinate-only detector was built (`tools/detect_overassembly.py`, no
library sequence needed) and run over all 31,078 families in the five
clusterable BEDs; 3,650 were testable, 69 flagged. **It does not support a
rate**, on three independent grounds:

1. **Calibration against a decoy null.** Applying the identical criterion at
   *half-integer* k′ (2.5 / 3.5 / 4.5), which no tandem array can produce,
   gives real 70 vs decoy 74.7 once matched — **ratio 0.94, p = 0.72**. Per
   tool: rm2 0.78, pantera 0.30, repet 1.16 (p = 0.17). There is no excess of
   integer-k structure anywhere.
2. **Measured sensitivity is ~0 on real cases.** An all-vs-all blastn of the
   rm2 consensus library — which was on disk the whole time — finds **68 of
   2,550 rm2 consensi carrying internal same-strand tandem structure**. The
   detector flags essentially none of them: sensitivity **0/8** on
   sequence-verified rm2 over-assemblies. `rnd-3_family-566` (821 bp) is a
   textbook 3 × 273 bp array by self-alignment (560 bp HSP at offset 266,
   92.5%; 287 bp at offset 522) and fails the detector's own length-
   commensurability condition. **So "rm2 shows no over-assemblies" is false** —
   the screen simply cannot see them.
3. **Tuning-dependence.** The phase histogram's bin count is arbitrary; three
   of the four confirmed families lose their flag at `n_bins=16`. And the
   enrichment level claimed as diagnostic (≥2.6) is reached by **6.1% of all
   testable families**, 194 of which are not flagged.

The cross-tool ratio that *does* confirm the four REPET entries is not
specific on its own either: across all 1,546 cluster members with
length-carrying mates, **243 hit some integer k in 2…5** — 103 rm2, 139 repet,
47 pantera. Confirmation needs both legs, and even then it is a screen.

**The honest statement:** four REPET entries in this assembly are confirmed
over-assemblies (chimeric, per the correction above). Nothing can be said about how common the phenomenon is per
tool, because the only sequence-free screen available has ~zero measured
sensitivity. *A per-tool rate requires the libraries.*

### The pipeline consequence, and its size

Judging copies of an over-assembled entry against its own consensus makes real
copies look 71% truncated. Pooled over the four confirmed families (10,229
hits), near-full-length under each rule, against a cross-tool ground truth of
**49.5%** taken from the independent rm2 consensus of the same element:

| rule | near-full-length |
|---|---|
| its own 764–771 bp consensus (old behaviour) | **0.74%** — ~67× too low |
| unit length, "both ends reached" test kept | 25.2% — wrong, mixes coordinate systems |
| **span-only against the unit length** | **58.2%** — closest to truth |

**Implemented:** when a member's modal consensus length is ≥1.8× the median of
its cluster-mates', stage 1 substitutes the mate length, records
`cons_len_source = deconvolved`, and switches that member to a **span-only**
near-full-length test (its `repeat_start`/`repeat_end` live in the collapsed
entry's coordinate system, so "reached both ends" is not a meaningful question
about one unit).

**Deliberately not divide-by-k**: k is not identified by coordinates, and of 18
flagged families with length-carrying mates only 4 had L/mate ≈ k. For the
other 14 the consensus was already *shorter* than their mates', so dividing
would have compounded the error up to 5×.

**And the rule abstains when the mates disagree with each other.** Inheritance
can otherwise inherit an over-assembly: rm2 `rnd-4_family-1503` (2,533 bp, a
sequence-verified dimer) has exactly one length-carrying mate, so a median over
mates would have been a single opinion, not a consensus. Abstentions are
recorded in `deconvolved.tsv` alongside the substitutions.

### Measured effect of turning it on — and it is not a REPET-only problem

Running the rule over the 12-cluster batch: **6 members deconvolved, 0
abstentions**. Near-full-length before → after, same copies, same gate:

| member | own consensus → mate reference | near-full-length | factor |
|---|---|---|---|
| `rm2:rnd-4_family-1870` (gap_aware) | 1,992 → 322 bp | 0.06% → **29.06%** | **494×** |
| `rm2:rnd-4_family-1870` (merge_always) | 1,992 → 322 bp | 0.06% → 32.68% | 509× |
| `repet:…G876-Map3` | 655 → 298 bp | 0.66% → 54.44% | 82× |
| `repet:…G1303-Map20` | 764 → 266 bp | 1.44% → 64.66% | 45× |
| `repet:…G1473-Map8` | 765 → 266 bp | 1.42% → 63.65% | 45× |
| `repet:…G1747-Map3` | 997 → 488 bp | 2.00% → 63.03% | 32× |
| `repet:…G2411-Map5` | 1,450 → 436 bp | 2.70% → 65.45% | 24× |

**The largest single case is rm2, not REPET** — `rnd-4_family-1870` at 1,992 bp
against cluster-mates at 322 bp, a ratio of 6.19. This matters because the
first pass of this investigation concluded that "rm2, pantera and fastltr show
no evidence at all", on the strength of a coordinate-only screen. That
conclusion was wrong twice over: refuted from the library side (68 of 2,550 rm2
consensi carry internal same-strand tandem structure) and now refuted from the
pipeline side, where the cross-tool ratio finds an rm2 entry with the biggest
effect in the batch. **Over-assembly is not a REPET quirk.**

The knock-on effect is larger than the seven rows above, because members that
*inherit* a length from cluster-mates inherit the corrected one: EDTA members
in the same clusters moved from 1,992 → 322 bp and 1,450 → 436 bp. Final
provenance across the batch: `bed16` 499,476 rows / `clustermate` 144,157 /
`deconvolved` 62,193 / `none` 141,287.

**Slide home:** stage 1.5 slide (`d9_deconvolution.png`) — DONE.
The strongest version of the slide is the 0.74% → 49.5% gap and the 494× row,
plus the negative result: the tempting sequence-free detector does not work,
and the thing that does work is *having the cluster* — which is the argument
for the whole pipeline.

---

## F5. Two tools cannot participate in family clustering (already in the deck)

`windowmasker`: **1 name for 4.36 M hits** — a single pseudo-family spanning
hundreds of Mb links to everything (inflated the graph to 9,461 edges and
lenient's largest cluster to 3,787 members before exclusion).
`ltrdenovo`: 3,104 names / 3,104 hits — per-copy ids, no family concept.
Both still count for coverage/support upstream.

**Slide home:** stage 0.1 (`d1_availability.png`) — already there.

---

## F6. Consensus length must be modal, not maximal

`rnd-1_family-286`: 31,387 hits imply a 269 bp consensus; **11 outlier hits**
imply 4,202 bp. Taking the maximum made near-full-length copies read **0**
where the true count is ~14,600. This is the `len_inconsistent` flag from the
report work biting the pipeline directly.

**Slide home:** stage 1.2 (boundaries) notes, alongside the median-not-max rule
for the over-extension guard — same lesson, two places.

---

## F7. Cluster members re-annotate the same loci 2.6x over — pooling must deduplicate

**Measured (stage 2, cluster 62, gap-aware):** the five members of the cluster
emit **90,313 copy rows** between them, but those rows describe only **35,194
distinct genomic loci** — a mean of **2.57 annotations per locus** (max 8).
Pairwise, the overlap is near-total: rm2 vs pantera **94.3%** of the smaller
member's copies, edta vs rm2 **91.8%**, edta vs pantera **88.4%**.

This is not a defect — it is the *definition* of a cluster; the members were
grouped precisely because they describe the same element. But it means a
pooled seed built by sampling copy rows would place the same sequence in the
alignment ~2.6 times on average, inflating apparent depth and weighting the
consensus toward whichever tool annotates most prolifically (rm2 contributes
59% of representatives, pantera 30%, edta 5%).

**Rule adopted:** collapse to one representative per locus before sampling.
Two annotations are the same locus when they overlap **reciprocally by >=50%
of the shorter** — the same test stage 0 uses for graph edges. Mere adjacency
is not enough: this element occurs in tandem arrays, and a touch-based rule
chained whole arrays into single "loci" (max members per locus 30 under
touching, **8** under reciprocal overlap). It is F1's failure mode a third
time, now at the locus level.

**The payoff — cross-tool support predicts copy completeness.** With loci
deduplicated, the number of *tools* supporting a locus can be counted, and the
near-full-length fraction rises monotonically with it:

| tools at locus | loci | fraction near-full-length |
|---|---|---|
| 1 | 5,126 | **0.4%** |
| 2 | 8,795 | 25.7% |
| 3 | 17,885 | 63.6% |
| 4 | 3,388 | **73.4%** |

A locus only one tool sees is almost never a complete copy; a locus all four
see usually is. This is the first direct measurement that multi-tool agreement
carries information about copy quality rather than merely about tool overlap,
and it is the quantitative justification for the whole clustering premise.

**Batch-level figures now on the slide** (12 clusters, gap-aware, not only cluster 62):
redundancy mean **2.17x** (range 1.26–3.04); pooled support/completeness **1 tool
108,884 loci at 5.9% · 2 tools 50,964 at 28.1% · 3 tools 43,214 at 54.0% · 4 tools
18,255 at 67.1%** — the monotonic rise holds across the whole batch.

**Slide home:** stage 2.1 slide (`d10_pooling.png`) — DONE.

---

## F8. At seed level the merge modes are indistinguishable — the decision falls to the secondary criteria

**Measured across 12 candidate clusters, both modes, 100 sampled copies each,
MAFFT** (rebuilt consensus lengths span 265–721 bp, so this is not one element
type):

| quantity | gap_aware − merge_always (median) | clusters where gap_aware is better |
|---|---|---|
| rebuilt consensus length | **0.0 bp** | 5 longer / 5 shorter / 2 equal |
| median seed depth | **+0.5** | 6 / 12 |
| deduplicated loci available | **+580** | **12 / 12** |
| bridged (multi-fragment) rows in the seed | **−4** | **11 / 12** |
| copies over 1.5× the modal consensus | **−6** | **8 / 12** (4 tied, 0 worse) |

**The test the plan expected to be decisive is not decisive.** Seed depth and
rebuilt consensus length — RESUME step 4.6's stated criteria — do not separate
the modes: the consensus differs by under a base pair in the median, and the
direction is not even consistent (4 vs 6). Whatever merge-always absorbs into
its copies, the alignment and the consensus caller discard again.

The modes separate cleanly on the secondary criteria, and gap-aware wins all
three: it yields **more loci in every single cluster** (12/12), puts fewer
bridged rows into the seed in **11 of 12**, and fewer copies over 1.5× the modal
consensus in **8 of 12** (tied in the other 4, never worse).

**Decision: `gap_aware` confirmed as the default** — not because it builds a
better consensus, but because it reaches the same consensus from cleaner, more
numerous, less chimeric evidence. **This supersedes the "still to be confirmed
at seed level" note in F3.**

*(Recomputed from `seed_mode_comparison.tsv` while building the deck: median consensus
difference 0.0 bp on a 5/5/2 split, median depth +0.5 better in 6 of 12; bridged rows
better in 11 of 12 and over-length copies in 8 of 12 (not 10/12 for either). The earlier
−0.5 bp / +1.0 / 7-of-12 figures did not reproduce; the conclusion is unchanged and
slightly stronger — depth and consensus length separate the modes even less.)*

**Slide home:** stage 1.1, replacing the provisional wording — DONE. The honest framing
is the interesting one: we expected depth to decide it, and depth had nothing
to say.

---

## F8b. Engine disagreement is NOT a triage signal — refuted at batch scale

**The hypothesis, from n=24:** the two packets where MAFFT and Refiner disagreed
most were both in the one cluster whose canonical classification path stopped at
`repeat:TE:ClassII`, suggesting `|len(mafft) − len(refiner)|` might flag the
clusters a curator should look at first.

**Tested properly on 412 packets from 207 clusters**, stratified across
classification depth and tool support rather than taken from the top of the
copy-number distribution (the original 12 were all top-2% by copy number and
9 of 12 were 4-tool — the stratum where tools agree by construction):

| path depth | packets | median &#124;Δ&#124; | median relative Δ |
|---|---|---|---|
| 3 (stops at Class I/II) | 66 | 72.0 bp | 9.8% |
| 4 (order) | 116 | 47.0 bp | 6.4% |
| 5 (superfamily) | 230 | 58.5 bp | 8.0% |

**Not supported.** Spearman(depth, |Δ|) = **+0.081, permutation p = 0.10**; on
the length-normalised measure **+0.060, p = 0.22**. The one marginal result
(Mann-Whitney on absolute Δ, p = 0.046) points the **opposite way** from the
hypothesis — shallow-path clusters disagree *less*, not more — and the effect is
non-monotonic across the three depths. Three tests, no consistent direction, no
effect that survives normalisation: this is a null result, and the original
observation was two packets in one cluster.

**What is true and worth keeping:** the engines disagree considerably more than
the first batch suggested — median **53.5 bp (7.7%)** of consensus length, with
exact agreement in only 2.7% of packets. On the original 12 clusters the median
was 4 bp. The difference is entirely composition: those 12 were short,
high-copy, well-covered families where any aligner converges. Engine choice
matters much more on the rest of the repeatome than the first sample implied —
which is an argument for running both engines, just not the triage argument F8b
originally made.

**Slide home:** keep the two-engine slide, drop the triage claim. The honest
framing — *we proposed a cheap curator-priority signal, tested it at 17x the
sample size, and it is not there* — is worth more to the deck than a
correlation that would not replicate.

---

## F9. `stk lint --genome` caught two identifier-integrity bugs that would have shipped silently

Both were invisible to every check that does not compare the alignment row
against the assembly, and both produce a seed that is *structurally* valid —
correct field set, correct widths, correct gap characters — while asserting
coordinates whose sequence is not what sits beside them.

**1. MAFFT silently reverse-complements rows.** `--adjustdirectionaccurately`
flips rows it judges backwards and marks them with an `_R_` prefix. Measured:
**77 of 100 rows** in cluster 540 gap-aware, 3 of 97 in cluster 62. The
sequence in a flipped row is the opposite strand from the one its Smitten
identifier names, so the strand must be flipped with it.

**2. EDTA leaves the strand uncalled.** `.` appears on **~0.3%** of EDTA hits
(4,189 of 1.39 M sampled) and on no other tool's. Extraction reverse-
complements only on `-`, so an uncalled strand is extracted *forward* — but
the identifier-retagging code treated `.` as reverse. The Smitten format has
no `.` anyway. Normalised where fragments are created, so extraction and
identifier can no longer disagree.

A third, milder issue in the same family: trimming the alignment removes
terminal bases, so a trimmed row no longer contains the whole interval its
identifier claims. Every row's coordinates are now walked in by the number of
non-gap bases actually removed from each end.

**Why it matters for the slides:** this is the concrete argument for running
`stk lint --genome` rather than `stk lint` alone, and for keeping lint as a
*visible soft gate*. All three bugs pass tier-1 validation cleanly. Dfam's
downstream TSD and extension algorithms retrieve flanking sequence **by
identifier** — a seed with subtly wrong coordinates is worse than one that
fails loudly, because it fails later, in someone else's pipeline.

**Also found:** the placeholder `au_string` "C. Goubert, Wheeler Lab" is a lint
**ERROR** (`au_format`) — abbreviated first names are rejected. A spelled-out
first name passes; an ORCID prefix downgrades the remaining INFO.

**Slide home:** stage 3 slide (`d11_seed_lint.png`) — DONE.

---

## F10. Consensi rebuilt from track copies reproduce the tools' own libraries

**Measured:** the stage-2 consensus, built from <=100 sampled genomic copies
with no access to any tool's library, aligned against RepeatModeler2's
independent consensus for the same family (GenomeArk
`RepeatModeler-v2.0.8/fasta/GCA_951799975.1-families.fa`):

| cluster | rebuilt | rm2 library | identity | aligned span |
|---|---|---|---|---|
| 62 | 265 bp | `rnd-1_family-286`, 269 bp | **98.1%** | 265/265 bp (full) |
| 540 | 363 bp | `rnd-2_family-317`, 378 bp | **95.2%** | 377 bp (full, minus strand) |

The rebuild is not merely *similar* to the library entry, it is co-terminal
with it: cluster 62's 265 bp maps to library positions 3–267 of 269.

### The classification the seed ships can contradict the sequence it ships

Cluster 62's rebuilt consensus carries a **terminal inverted repeat** — 14 bp
with one mismatch, the first 10 bp (`CTTTAAAGGG`) perfect — and is 265 bp with
no coding capacity. That is a MITE. The cluster's member-weighted majority
path is `repeat:TE:ClassI:LINE:Jockey`, so the emitted seed carries
**`#=GF TP LINE/I-Jockey`**: a Class I retrotransposon label on a Class II
element.

This is the F2 caveat made concrete and consequential. The tool labels are
wrong here, the pipeline faithfully propagates them into a Dfam-format field,
and nothing in the lint chain can object — `tp_unknown` only checks that the
value exists in Dfam's vocabulary, not that it fits the sequence.

**Consequences worth acting on:**
1. The curator queue must show the rebuilt consensus and its structural
   features *next to* the inherited classification, and treat the
   classification as a proposal, not a datum. The verdict schema already
   carries a classification field; this is why.
2. A cheap automated cross-check is available: a TIR at both ends of a short
   consensus contradicts any Class I `TP`. Class I vs Class II is the one
   split where structure is decisive enough to flag automatically, and TE-Aid
   already computes what is needed.
3. Until then, `#=GF TP` should be read as "what the tools said", and that
   should be stated in the packet rather than implied.

**Slide home:** validation slide (`d12_validation.png`) — DONE. The counterweight: the rebuild works,
and the label attached to it may still be wrong.

**Every seed row checked against the assembly.** Across the 12-cluster batch
(48 seeds, both engines, both merge modes), **5,005 of 5,005** alignment rows
reproduce the genome exactly at the coordinates their Smitten identifier
claims, and `stk lint --genome` reports **0 ERRORs and 0 coordinate warnings**.
The check is done independently of lint because lint's own coordinate checks
are advisory (F9).

**Against RepeatModeler's own seed alignment.** GenomeArk also publishes
rm2's `.stk` seeds for this assembly, so the comparison can be made at the
*seed* level rather than only at the consensus level. Its record for
`rnd-1_family-286` carries **100 sequence rows, median span 265 bp** — the same
length this pipeline rebuilds. **All 100 of those loci are present in our
locus set** (>=50% reciprocal overlap), which contains **35,205**. So the
track-derived locus discovery is a strict superset of what RepeatModeler
sampled for its own seed, and the disagreement between the two seeds is a
matter of *which* copies were chosen, not of which copies exist.

Worth noting for the deck: rm2's own seed record would **not** pass
`stk lint`. It carries `DE`, `TP`, `CC`, `BM` and `SQ` but no `AU` and no `OC`
(both required, both ERRORs), and its identifiers omit the assembly accession
(`OX637605.1:2145865-2146095_+` rather than Smitten). The gap between "what a
tool emits" and "what Dfam accepts" is exactly the gap this pipeline closes.

**Slide home:** the validation slide — this is the evidence that rebuilding
from copies works at all, paired with the seed-level comparison above.

---

## F11. An LTR-only tool in a cluster is an independent quality marker — the structural finders are worth their own lane

**Prompted by:** the observation that libraries exist for the LTR-specific
finders too, and that LTR finders may bear low-hanging fruit in combination
with — or independently of — the general-purpose repeat finders.

**fastltr is extraordinarily specific.** Across the 1,206 strict clusters it
joins only **117** (9.7%, versus rm2's 92.9% and edta's 89.1%), and of those
**113 are `LTR` by majority order — 96.6%**. Three are TIR, one Helitron. It
never appears alone: edta co-occurs in 116/117 clusters and rm2 in 112/117.
It is not a general-purpose participant; it is a *witness*.

**Its presence marks better clusters, and not merely bigger ones.** Among all
LTR clusters, those containing fastltr pass the candidate gates at **89% vs
56%**:

| gate | with fastltr | without |
|---|---|---|
| G1 full-length | 0.94 | 0.69 |
| G2 depth | 0.99 | 0.78 |
| G4 tandem | 1.00 | 0.93 |
| **candidate (all gates)** | **0.89** | **0.56** |

fastltr necessarily adds a tool, so the obvious confound is cluster size.
Controlling for it by **matching on exactly which other tools are present**,
the effect survives:

| other tools present | candidate rate with fastltr | without | n (with/without) |
|---|---|---|---|
| edta+pantera+repet+rm2 | **0.97** | 0.70 | 33 / 30 |
| edta+repet+rm2 | **0.97** | 0.69 | 29 / 52 |
| edta+rm2 | **0.83** | 0.48 | 29 / 115 |

In the top stratum the fastltr-supported clusters are *smaller* (median pooled
193–397 kb vs 821 kb) and still pass more often, so this is not a size effect.

**Reading.** fastltr finds elements by *structure* — terminal repeats and
target-site duplications — not by homology to a library. When a
structure-based caller and several homology-based callers converge on the same
locus, the agreement is across **methods**, not merely across implementations
of the same method. That is a qualitatively stronger form of support than one
more homology tool agreeing, and it is measurable here.

**A caveat that runs the other way:** perfect order coherence is *lower* among
fastltr clusters (75% vs 92% at coherence = 1.0), purely because adding a
fourth or fifth voter adds an opportunity to disagree. Order coherence is
comparative and member-weighted (F2); it should not be read as fastltr making
clusters messier.

**`ltrdenovo` is a second structural channel, currently unused.** 3,104 hits,
one name per hit, so it cannot cluster (F5) — but every hit is classified LTR
(RLG 1,133 / RLC 1,109 / Unknown 861), median 1,310 bp, 7.1 Mb total. As
per-copy structural predictions these are exactly the kind of orthogonal
evidence that could corroborate an LTR cluster without ever joining it.

**Proposed pipeline consequence (not yet implemented):** treat structural
LTR evidence as a *support channel* rather than a clustering participant —
a per-cluster flag `structural_ltr_support` counting overlapping fastltr
members and ltrdenovo copies. On this assembly it would promote ~100 LTR
clusters that are already candidates, and — more usefully — mark which of the
119 fastltr-free LTR clusters lack structural corroboration and should be
looked at harder.

**Open question for the group:** LTR elements are the easiest class to
validate by eye (LTRs, TSD, PBS/PPT, and a rebuilt consensus that should be
co-terminal with the LTRs). If "low-hanging fruit" means *cheapest to curate
confidently*, the fastltr-supported LTR clusters may be a better first batch
than the highest-copy-number families, which are dominated by short
non-autonomous elements.

**Slide home:** a new stage-0 slide on tool scope — pairs naturally with
`d1_availability.png` and F5 (which tools *cannot* cluster) by adding which
tools *should not be expected to*, and what they are good for instead.

---

## F12. A seed with no consensus lints clean — the batch found 184 of them

**Measured (first 200-cluster batch, 414 packets):** **174 of 414 Refiner
packets and 10 MAFFT packets** produced a rebuilt consensus of length **zero**,
and every one of them emitted a `.stk` file that passes `stk lint --genome`
with no ERROR and no WARN. An all-`.` RF line is exactly as wide as its
all-gap sequence rows, so nothing in the format is violated. Structurally
valid, semantically empty — the same shape of failure as F9, one level up.

**It was not the aligner's fault.** Refiner had called a consensus in
**174/174** of those packets (median 4,217 bp, max 15,925 bp) and every one of
its sequence identifiers parsed. The defect was in this pipeline's occupancy
rule: a column counted as part of the consensus only if ≥50% of **all** rows
were non-gap there. For a 4.2 kb element whose copies are fragments of
different regions, no column ever reaches 50% of all rows, so the alignment has
no match columns and the consensus is empty.

**Fix — occupancy is span-normalised.** A column's denominator counts only the
rows that *reach* it: a copy ending at position 900 says nothing about position
3000 and should not be counted against it. Paired with an absolute depth floor
(≥3 rows), because the two catch opposite failures — span-normalisation alone
keeps a ragged shoulder that only two rows reach (2/2 reads as fully occupied),
and the depth floor alone is what discarded the 174.

**Validated against an independent reference.** Both clusters whose consensus
had already been checked against rm2's own library move *toward* it:

| cluster | before | after | rm2 library | identity before → after |
|---|---|---|---|---|
| 1103 | 265 bp | **270 bp** | 269 bp | 98.87% → **98.89%**, now full length |
| 1048 | 363 bp | **380 bp** | 378 bp | 95.49% → **98.68%** |

**After the fix: 0 of 414 packets degenerate**, and 88,368 alignment rows across
both engines verify byte-exact against the assembly.

**Why the first batch could never have found this.** The original twelve
clusters were all short, high-copy, well-covered families in the top 2% of the
copy-number distribution — every row spanned the whole element, so a 50%
all-rows threshold was always attainable. The defect only exists for long
elements with fragmentary copies, which the stratified batch was built to
include. **A sample chosen for being easy to build cannot test the builder.**

A packet with no match column is now marked `degenerate` and raises a
`degenerate_alignment` ERROR in `lint_triage.tsv`, because `stk lint` cannot
see it and never will.

**Slide home:** pairs with F9 as the second half of "what a format checker
cannot see" — F9 is coordinates that are wrong, F12 is a consensus that is
absent.

---

## F13. A cluster's members do not always describe one length — the anti-tandem guard deleted the full-length evidence

**Prompted by:** CG's review of the curator demo — *the rebuilt consensus is
systematically no better than the best member, and always lands near the
shortest one.*

**Measured (cluster 1183, an LTR/Gypsy family).** Its eight members carry two
quite different lengths: three at 414–423 bp and five at 6,523–7,400 bp. The
pipeline reduced that to one modal length, picked **423 bp**, and the F4
anti-tandem cap (1.5 × modal) then excluded **all 83 loci above 635 bp** — 52
REPET, 12 rm2, 9 fastltr, 8 edta, i.e. every full-length copy. Sampling saw
only the short side (max sampled span 618 bp) and the seed came out at
**453 bp for a 7.4 kb element**.

Two compounding causes, both worth stating:

1. **The mode was taken over copy ROWS, not over deduplicated loci.** A member
   with many redundant annotations decides the cluster's length. Here 423 bp
   wins on 326 rows against REPET's 189 at 7,400 — and only because
   `edta:TE_00000385` contributes no length of its own. Over deduplicated loci
   the mode is 7,400.
2. **The F4 guard is right in general and wrong here.** It was tuned on a
   265 bp element whose REPET entry was a chimera (F4); on a bimodal cluster
   the same rule deletes the real element instead of the artifact.

**The fix, and why it does not reopen F4.** Member lengths are grouped into
**modes** (split where consecutive lengths differ by >2×), and a mode counts
as real only when **≥2 distinct tools measured it**. That is the same
cross-tool principle the deconvolution guard uses, and it separates the two
cases cleanly:

| cluster | short mode | long mode | outcome |
|---|---|---|---|
| 1183 | 418 bp, 2 tools | 7,363 bp, **4 tools** | both real → two models |
| 62 (F4) | 264 bp, 2 tools | 764 bp, **REPET alone** | long mode rejected → F4 holds |

Each locus is then capped against **the mode it belongs to**, nearest in log
space — so the boundary sits at the geometric mean (~1,765 bp for 1183). A
tandem dimer of the solo LTR at 850 bp is still capped; a truncated full
element at 1.8 kb is judged against 7,363 and kept.

**Only MEASURED lengths vote.** A member whose `cons_len_source` is
`clustermate` or `deconvolved` is repeating a cluster-mate's opinion. Counting
its tool inflates the mode's apparent support — and since `#=GF CC` names the
supporting tools, it would ship a false provenance claim. On 1183, EDTA
contributes 423 and 7,400 bp while **EDTA's own library entries are 244 and
272 bp** and `edta.bed` carries `NA` consensus coordinates: EDTA never measured
either length.

**Result on 1183:**

| model | MAFFT | Refiner | near-full-length | lint |
|---|---|---|---|---|
| `model0_418bp` (solo LTR) | 455 bp | 426 bp | 44% | clean |
| `model1_7363bp` (full element) | **7,389 bp** | **7,372 bp** | **97%** | clean |

**Scope:** 12 of 207 clusters have two real length modes, **9 of them LTR**.
Only 1183, 1201 and 1162 had actually landed on the short mode; the other nine
happened to pick the long one already, so the damage was invisible in aggregate
and only showed up because CG looked at the flagship packets by eye.

**Slide home:** stage 2, next to F4 — the same guard, the case where it helps
and the case where it hurts, and the cross-tool test that tells them apart.

---

## F14. Cluster 1183 is a Ty3/Gypsy element, decomposed from sequence — and one member does not belong

Worth recording because it is the first time the pipeline's output has been
checked against element *structure* rather than against another tool's
consensus, and because it settles a classification contradiction.

**The decomposition is exact.** rm2's 423 bp `ltr-1_family-35` (its own
`Type=LTR` model) matches **both ends and only the ends** of the long consensi.
REPET's 7,400 bp entry carries LTRs at 12–437 and 6,965–7,390 — **426 bp and
100% identical to each other** — both beginning `TG` and ending `CA`. rm2's
6,523 bp `Type=INT` model fills 438–6,964 at 99.77%. The arithmetic closes:
**6,523 + 2 × 423 = 7,369**.

Supporting evidence: **57 of 60** full-length genomic loci carry a 4 bp target
site duplication (and none carry 5 or 6), independently reproducing FastLTR's
`tsdl4_tsdc0.94` header; `getorf` finds exactly one ORF, 1,898 aa with no
internal stops, lying entirely inside the internal region and containing a GAG
CCHC zinc knuckle and the RT catalytic `YLDD`; FastLTR annotates
`domains=PROT|Ty3_gypsy RT|Ty3_gypsy RH|Ty3_gypsy INT|Ty3_gypsy`.

**So `majority_path = repeat:TE:ClassI:LTR:Gypsy` is correct** and the EDTA
`MITE/DTA` label is not a counter-argument — see F15.

**RETRACTED — `pantera:Gypsy_9` is a genuine member.** An earlier draft of this
entry claimed it was mis-clustered, on the strength of a blastn showing zero
homology to the other members. That was an artefact of the wrong input file.

Tested the decisive way instead — take the GENOMIC sequence at a locus pantera
labelled `Gypsy_9` and ask what it is. At
`OX637597.1:39,330,768-39,338,134` (7,367 bp):

| against | identity | span |
|---|---|---|
| fastltr `CONS_4-7367` | **98.96%** | 7,368 bp, full length |
| repet `G3000-Map20` | 98.98% | 7,377 bp |
| rm2 `ltr-1_family-36` (INT) | 99.05% | 6,530 bp |
| rm2 `ltr-1_family-35` (LTR) | 98.58% | 423 bp |
| `Gypsy_9` in the pantera library on disk | **no hits** | — |

**The library file does not correspond to the BED.**
`GCA_951799975.1.fGobNig-pantera-pass.fa` is a different pantera run: 196
sequences against 14,325 family names in `pantera.bed`, only 135 names in
common, and where a name occurs in both the sequences are unrelated.
`Gypsy_9-fGobNig` exists in both and they are different elements.

**Two lessons worth keeping.** First, a per-tool sequence analysis is only as
good as the provenance of the library — F4's per-tool rates need the library
that actually produced the BED, and for pantera we do not have it. Second, the
right test for "is this member really this family" is against the GENOME at the
member's own annotated loci, not against a library entry that shares its name.

**Slide home:** the validation slide — this is what "rebuilt from copies"
delivers when checked against structure, plus the honest caveat that clustering
is coordinate-based.

---

## F15. EDTA's library header names a different sequence from the one it annotates

**Terminology, per CG:** *MITE is not a classification label — it is a
qualifier for a non-autonomous DNA/TIR element.* So `MITE/DTA` reads as
"non-autonomous DTA (hAT)", and the class is DNA/TIR-hAT. The entries below are
therefore correctly *qualified*; the problem is not the qualifier.

**Measured.** EDTA's `MITE/` qualifier tracks length closely: 1,956 such
entries, median 292 bp, **maximum 599 bp** — consistent with non-autonomous
elements, which are short by definition. Cluster 1183's EDTA members are 244 bp
and 272 bp in the library, so calling them non-autonomous TIR elements is
defensible *for those sequences*. What is not defensible is that EDTA reuses
the same family name for ~7.4 kb structural LTR loci in its annotation: the
library entry and the annotated loci are not the same thing.

**The header disagrees with EDTA's own BED annotation on 27.3% of families.**
Across the 207 built clusters, EDTA member classes contradict the cluster's
majority order for **46 of 203 members (22.7%)**, in 39 of 164 clusters;
restricting to clusters with a resolved (non-`TE`) order, **22 of 179 (12.3%)**,
dominated by LTR→TIR.

**Consequence:** the pipeline already takes class from the BED rather than the
library header, so nothing is broken — but this is the measurement that says
*never* switch to the header, and the residual ~12.3% genuine BED-class
conflict is a triage list rather than a bug.

**Slide home:** a notes item beside F5 (which tools can participate) — this is
about which *fields* of a participating tool can be trusted.

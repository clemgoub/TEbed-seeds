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

## F4. A tool consensus can be an over-assembly of tandem units — confirmed for four REPET entries, and NOT generalisable

**Prompted by:** cluster 62's two REPET families having suspiciously clean
consensus lengths of **764 and 765 bp**, against 264–269 bp from rm2, pantera
and edta for the same element, and the hypothesis that this was a
palindromic/MITE structure worth tolerating.

### Confirmed: the 765 bp entries are ~3 tandem units of a ~255 bp element

- **Palindromic structure: refuted.** Adjacent rm2 hits within 700 bp are
  **75% same-strand**; a palindrome would show an *excess* of opposite-strand
  pairs, not a deficit.
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
over-assemblies. Nothing can be said about how common the phenomenon is per
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

**Slide home:** a stage-1.5 slide, "when the library consensus is the problem".
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

**Slide home:** a new stage-2 slide on pooling — the redundancy figure and the
support/completeness table belong together.

---

## F8. At seed level the merge modes are indistinguishable — the decision falls to the secondary criteria

**Measured across 12 candidate clusters, both modes, 100 sampled copies each,
MAFFT** (rebuilt consensus lengths span 265–721 bp, so this is not one element
type):

| quantity | gap_aware − merge_always (median) | clusters where gap_aware is better |
|---|---|---|
| rebuilt consensus length | **−0.5 bp** | 4 better / 6 worse / 2 equal |
| median seed depth | **+1.0** | 7 / 12 |
| deduplicated loci available | **+580** | **12 / 12** |
| bridged (multi-fragment) rows in the seed | **−4** | **10 / 12** |
| copies over 1.5× the modal consensus | **−9.5** | **10 / 12** |

**The test the plan expected to be decisive is not decisive.** Seed depth and
rebuilt consensus length — RESUME step 4.6's stated criteria — do not separate
the modes: the consensus differs by under a base pair in the median, and the
direction is not even consistent (4 vs 6). Whatever merge-always absorbs into
its copies, the alignment and the consensus caller discard again.

The modes separate cleanly on the secondary criteria, and gap-aware wins all
three: it yields **more loci in every single cluster**, and in 10 of 12 it puts
fewer bridged rows and fewer over-length copies into the seed.

**Decision: `gap_aware` confirmed as the default** — not because it builds a
better consensus, but because it reaches the same consensus from cleaner, more
numerous, less chimeric evidence. **This supersedes the "still to be confirmed
at seed level" note in F3.**

**Slide home:** stage 1.1, replacing the provisional wording. The honest framing
is the interesting one: we expected depth to decide it, and depth had nothing
to say.

---

## F8b. Engine disagreement is a triage signal, not just noise

**Measured:** the same 24 packets built twice from the *same* sampled copies,
once with MAFFT and once with Dfam's Refiner. Consensus lengths agree to a
median of **4 bp** (identical in 5 of 24, within 5 bp in 15 of 24). MAFFT is
the longer of the two in 13 packets, Refiner in 6.

But the tail is not small: the two worst disagreements are **140 bp and 70 bp**,
and both are the same cluster — the one whose member-weighted classification
path stops at `repeat:TE:ClassII`, i.e. the tools could not agree what kind of
Class II element it is. The clusters where two independent aligners disagree
about how long the element is are the clusters where the tools also disagree
about what it is.

**Implication:** `|len(consensus_mafft) − len(consensus_refiner)|` is cheap
(both engines already run) and looks like a useful **queue-priority signal** —
a packet where the engines disagree is a packet where a curator's time is
well spent. It costs nothing to record and is worth testing against the first
batch of curator verdicts.

**Slide home:** stage 2/3, as the argument for running both engines rather than
picking one.

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

**Slide home:** stage 3 (seed + lint) — "what a format checker cannot see".

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

Cluster 62's consensus also carries a **12 bp terminal inverted repeat**
(self-alignment: positions 1–12 vs 267–256, minus strand, 91.7% identity),
i.e. it is a TIR element / MITE — although the cluster's member-weighted
majority order is `LINE`. That disagreement between what the tool labels say
and what the rebuilt sequence shows is exactly the F2 caveat, and it is an
argument for putting the rebuilt consensus, not the inherited label, in front
of the curator.

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

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

## F4. A tool consensus can be an over-assembly of tandem units — measured instance

**Prompted by:** the observation that cluster 62's two REPET families have
suspiciously clean consensus lengths of **764 and 765 bp**, and the hypothesis
that this reflects a palindromic/MITE structure worth tolerating.

**Palindromic structure: NOT supported.** Among rm2 hits of the cluster's
269 bp family, consecutive hits within 700 bp are **75% SAME strand** — only
**25% opposite**, i.e. *below* the 50% expected by chance. A palindromic
element whose library consensus captured one arm would produce an excess of
opposite-strand adjacent pairs. There is none.

**What the 765 bp actually looks like:**

- REPET copies of the 765 bp consensus have median length **223 bp** — they
  cover only **29% of their own consensus** (p75 = 266 bp = 35%).
- Their match coordinates along the consensus cluster into **three blocks**
  (start quantiles 2 / 264 / 539; end quantiles 262 / 511 / 765), i.e. three
  consecutive ~250 bp segments rather than a uniform spread.
- 765 / 3 ≈ 255 bp, and the independent rm2 and pantera consensi for the same
  cluster are **269 and 264 bp**.
- rm2 hits form **2,892 same-strand tandem runs** (2,727 of them exactly two
  units), median run span **623 bp**.

**Reading:** the true repeated unit is ~265 bp, and REPET's 765 bp entry is
most consistent with **~3 tandem units collapsed into one library consensus** —
the *consensus-level* analogue of F1's runaway merge chain. Both REPET families
in the cluster show the same length, and they partition loci almost perfectly
(0.2% overlap of 3,073 vs 3,280 copies), which is the classic signature of
redundant library entries dividing copies between them rather than
double-annotating.

**The implication inverts the original instinct:** the clean 765 bp pattern is
meaningful, but it should **not** be tolerated as the element length. Treating
it as ground truth makes every genuine copy look 71% truncated and would sink
them at any near-full-length gate. This is a concrete instance of the
"tool consensus may itself be wrong" caveat that has governed this project from
the start — now with a measurable signature: *median copy length ≈ consensus
length ÷ small integer, plus block-structured match coordinates.*

**Definitive test needs sequence** (not yet on hand): self-align the 765 bp
consensus — three tandem blocks would appear as off-diagonal repeats — and
align it against the rm2/pantera 265 bp consensi. Requires the REPET library
FASTA. **This is the parallel-agent task briefed in `RESUME.md`.**

**Slide home:** a new stage-1.5 slide, "when the library consensus is the
problem" — it motivates rebuilding consensi from copies, which is the whole
premise of the pipeline.

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

**Measured (stage 2, clusters 62 and 540, 100 sampled copies each):**

| cluster | mode | loci | rebuilt consensus | match columns | median depth | bridged rows | copies > 1.5x modal |
|---|---|---|---|---|---|---|---|
| 62 | merge_always | 34,738 | 265 bp | 265 | 73 | 2 | 588 |
| 62 | gap_aware | 35,194 | **265 bp** | 265 | 73 | **0** | **571** |
| 540 | merge_always | 11,486 | 364 bp | 364 | 66 | 5 | 112 |
| 540 | gap_aware | 11,978 | **363 bp** | 363 | 68 | **1** | **90** |

**The decisive test was expected to be seed depth and rebuilt consensus length
(RESUME step 4.6). It is not decisive: those two numbers are identical.**
Rebuilt consensus differs by 0 bp (cluster 62) and 1 bp (cluster 540); median
depth differs by 0 and 2. Whatever merge-always absorbs into its copies, the
alignment and consensus caller discard again.

The modes separate only on the secondary criteria, and there gap-aware wins on
every one at no measurable cost: **more loci** (+456, +492), **fewer bridged
(multi-fragment) rows** in the seed (0 vs 2, 1 vs 5), and **fewer over-length
copies** to exclude (571 vs 588, 90 vs 112).

**Decision: `gap_aware` confirmed as default** — not because it builds a better
consensus, but because it reaches the same consensus from cleaner, more
numerous, less chimeric evidence. Caveat: n=2 clusters; a wider run is under
way.

**This supersedes the "still to be confirmed at seed level" note in F3.**

**Slide home:** stage 1.1, replacing the provisional wording.

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

**Slide home:** the validation slide — this is the evidence that rebuilding
from copies works at all, and it is the natural place to introduce the
comparison against RepeatModeler's own seeds (its `.stk` for the same
assembly is public and carries 100 rows for `rnd-1_family-286`).

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

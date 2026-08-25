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

**Still to be confirmed at seed level (step 3):** the decisive comparison is
seed depth and rebuilt consensus length, not copy counts. **Update this entry
and the slide once the first seeds are built both ways.**

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

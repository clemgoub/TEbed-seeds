# Handoff — LTRDeNovo tested as a support channel: no F11-like effect (2026-08-26)

**Prompted by CG:** F11 shows fastltr's presence in a cluster lifts the
candidate pass rate 56% → 89%. *"What about LTRDeNovo? Not tested?"* Correct —
it was not. This file records the test and its result.

## Why it needed a different test

Per F5, `ltrdenovo` emits **3,104 names for 3,104 hits** — per-copy
identifiers, no family concept — so it is excluded from family clustering and
can never be a cluster *member* the way fastltr is. F11's test (member
presence vs candidate rate) is therefore not available for it.

Tested instead by **coordinate corroboration**: for each LTR cluster, what
fraction of its members' annotated intervals are ≥50% covered by some
ltrdenovo call. A cluster counts as corroborated above a coverage threshold.

## Result — no effect at any threshold

Strict linkage, 327 LTR-majority clusters. Reference row is F11 reproduced
exactly on the same table.

| test | n with | rate with | n without | rate without | delta |
|---|---|---|---|---|---|
| **fastltr membership (F11 reference)** | 113 | **0.894** | 214 | 0.556 | **+0.338** |
| ltrdenovo corroboration ≥2% | 86 | 0.663 | 241 | 0.676 | −0.014 |
| ltrdenovo ≥5% | 68 | 0.677 | 259 | 0.672 | +0.005 |
| ltrdenovo ≥10% | 58 | 0.655 | 269 | 0.677 | −0.021 |
| ltrdenovo ≥25% | 24 | 0.708 | 303 | 0.670 | +0.038 |
| ltrdenovo ≥50% | 2 | 0.500 | 325 | 0.674 | −0.174 |

Controlling for fastltr (the obvious confound), the effect is **negative** in
both strata: without fastltr 0.364 (n=22) vs 0.578 (n=192); with fastltr
0.833 (n=36) vs 0.922 (n=77). Gate-level detail at ≥10% shows no gain either
(G1 0.76 vs 0.78, G2 0.81 vs 0.86, G4 1.00 vs 0.94).

The two tools mark largely **different** clusters — of 327 LTR clusters,
fastltr joins 113, ltrdenovo corroborates 58, and only 36 overlap. So this is
not redundancy between two structural callers.

## Reading, and what would make it a fair test

The plausible interpretation is that **family-level participation**, not
coverage, is what F11 measures: a tool that groups copies into a family and
then agrees with other tools on that family is making a much stronger claim
than a per-copy caller whose interval happens to overlap. A per-copy caller
adds coverage without adding concordance.

This is a **first pass and should not be written up as F16 as it stands**:

1. The ≥50%-coverage / ≥X%-of-members criterion is invented here, not
   validated. ltrdenovo's median call is 1,310 bp against fastltr's
   full-element models, so partial coverage may be the norm rather than a
   signal of weak support.
2. 3,104 hits vs fastltr's ~45k makes ltrdenovo a much sparser witness;
   low power is not excluded as an explanation for the null.
3. **The decisive test is available**: group ltrdenovo's per-copy calls into
   families first (cluster its call sequences, or use its own internal
   grouping if the tool exposes one), then re-run F11's membership test.
   If a grouped ltrdenovo reproduces the fastltr effect, the conclusion is
   "structural detection is the signal"; if it does not, the conclusion is
   "fastltr specifically is the good witness" — those are different claims
   with different consequences for the support-lane design.

Data: `ltrdenovo_support_test.tsv` (the table above) and
`ltrdenovo_per_cluster.tsv` (per-cluster corroboration fractions and gates)
were produced in the deck session and can be regenerated from
`candidates_strict.tsv` + `inputs/ltrdenovo.bed` + the per-cluster `loci.tsv`.

**Deck status:** slide 18 now names fastltr explicitly and states the
ltrdenovo null with the "participation, not coverage" reading; the speaker
notes carry the numbers for the question CG asked, flagged as a first pass.

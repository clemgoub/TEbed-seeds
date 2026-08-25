# Next instructions — from CG's deck review (2026-08-25)

Read with `RESUME.md` (state) and `PIPELINE_FINDINGS.md` (measured results, and
the file the slides are rebuilt from). Branch: `stage2-seeds`.

You are idle waiting on two inputs. **Sections A and B need neither** — do them
now. Sections C and D unblock when the inputs arrive.

---

## Inputs still outstanding

| input | needed for | status |
|---|---|---|
| REPET consensus library FASTA | §C — sequence-level confirmation of over-assembly, per-tool rate | CG will provide |
| ORCID + spelled-out author name | §D — the clean re-run; `au_string` is currently a lint ERROR | CG approved the format |

`config/pipeline.yaml: au_string` is `"Clement Goubert"`. CG's decision:
**spelled-out first name plus ORCID prefix**, i.e.
`#=GF AU  ORCID:0000-0000-0000-0000 Clement Goubert` once he supplies the id.
Do not invent one. Until then the `orcid_missing` INFO is expected and fine.

---

## A. Scale the batch — highest value, no new inputs (do this first)

The engine-disagreement signal (F8b) is the finding CG most wants pursued, and
right now it rests on **6 discordant packets out of 24**. That is an anecdote.
This genome can supply far more without any new data:

- **763 strict candidates exist; 12 are built (1.6%).**
  Orders available: TIR 261 · LTR 220 · LINE 190 · TE (ClassII-unresolved) 34 ·
  Helitron 26 · PLE 15 · DIRS 9 · SINE 8.
- Measured cost: **17 s per cluster per engine pair** → a 200-cluster batch is
  **~2 h single-threaded**, less with `--threads`.

**Task A1 — stratified 200-cluster batch.** Sample across `majority_order` and
across `n_tools` (2/3/4/5), not just the top of `pooled_full_len` — the current
12 are all high-copy and 4-tool, which is exactly the stratum where tools
agree. Include the 34 `TE`-order clusters deliberately: those are where the
classification path stops shallow, and F8b's hypothesis is that they are also
where the engines disagree.

**Task A2 — test F8b properly.** Per packet record `|len(mafft) −
len(refiner)|` and the depth of `canonical_path` (3 = stops at ClassI/II,
5 = superfamily). Then test the association: is |Δ| higher where the path is
shallower? Report n, effect size, and a p-value from a rank test. Two honest
outcomes are both publishable-internally: it holds, or it was one cluster.

**Task A3 — write the result into `PIPELINE_FINDINGS.md` F8b**, replacing the
n=24 numbers, and say explicitly what one genome still cannot settle
(cross-taxon generalisation; calibration against real curator effort).

**Task A4 — extend F7's support/completeness table to the new batch.** It
currently pools 12 clusters (1 tool 5.9% → 4 tools 67.1% near-full-length).
Confirming that monotonic rise at 200 clusters makes it a headline result for
the VGP report's technical level, not just a pipeline internal.

---

## B. Didactic debt — CG's main criticism of the deck

CG's repeated note across four slides: *"I understand the intuition, not the
implementation"*, and *"real examples would be good here"*. The slides are
persuasive to someone who already knows the code. Fix that with real data.

**Task B1 — locus-dedup worked example (highest priority).** CG on the pooling
slide: *"it's not 100% explained how the representative is built."* Emit, for
ONE real locus in a built cluster, a small table: every member row (tool,
family, chrom, start, end, strand, n_fragments), which row became the
representative, the reciprocal-overlap value against each other row, and the
resulting support count. Save as
`runs/<asm>/examples/locus_dedup_example.tsv` plus a to-scale figure. State in
the docs why it is NOT the longest row (same reason boundaries do not come from
the widest member claim — it would let the most over-extended tool decide).

**Task B2 — one real edge for the graph slide.** For a single concrete edge:
both family names, both footprints, the joint bp, both weight channels
(`joint/min(footprint)` and copy co-annotation fraction), and one locus where
they overlap. → `examples/edge_example.tsv`.

**Task B3 — one real deconvolution case.** For `rm2:rnd-4_family-1870` (the
494× case): its own modal consensus length, each cluster-mate's, the ratio, the
substituted length, and near-full-length before/after. →
`examples/deconvolution_example.tsv`.

**Task B4 — split the stage-1.5 slide's material.** CG: *"I only understand
clearly the first plot."* Produce the numbers for a two-slide version — one
slide establishing that one family's consensus is an over-assembly, a second
generalising to the batch effect. Do not redraw the slides; produce the tables
and say what each panel should show.

---

## C. When the REPET libraries arrive

Context: four REPET entries are **confirmed** over-assemblies from coordinates
alone, but the sequence-free detector has **0/8 sensitivity** on
sequence-verified rm2 cases, so no per-tool rate can be claimed. See F4.

**CORRECTED 2026-08-25 — the mechanism is a CHIMERA, not tandem units.** CG
checked the locus in the genome browser and the tandem reading was wrong; see
F4's correction. Co-annotation identity per consensus block shows
`G1473-Map8` = rm2 `rnd-1_family-608` (496 bp) + rm2 `rnd-1_family-286`
(269 bp), and 496 + 269 = 765. `G1303-Map20` is the same pair fused in the
opposite order. Worked example, matching the browser:
`OX637613.1:25,978,734-25,979,682`, table at
`runs/<asm>/examples/overassembly_locus_example.tsv`.

**Task C1 — align each suspect consensus against its CLUSTER-MATE consensi**
(not against itself). The chimera predicts two non-overlapping, near-full-length
HSPs: `rnd-1_family-608` matching positions ~1-496 and `rnd-1_family-286`
matching ~497-765, each at high identity over most of the mate's length. A
self-alignment predicting internal same-strand repeats at ±255 is the TANDEM
hypothesis and should now be expected to FAIL — run it anyway as the control,
because if it succeeds the correction is wrong.

**Task C2 — align the long consensus against the short cluster-mates**
(rm2 `rnd-1_family-286` 269 bp, pantera `Unknown_572-fGobNig` 264 bp).
Over-assembly predicts the short one matching the long one ~3 times, at the
block boundaries (starts 2/264/539, ends 262/511/765).

**Task C3 — per-tool over-assembly rate, from sequence.** Two distinct
phenomena, so report them separately: (a) CHIMERAS — a consensus that aligns
end-to-end to two or more *different* shorter consensi covering disjoint parts
of its length; (b) TANDEM over-assemblies — internal same-strand repeat
structure. This is the number F4 currently cannot state. Note that the
coordinate-only phase detector conflates the two (496 ≈ 2 x 250 mimics a k=3
array), so sequence is required to separate them.

**Task C4 — validate the 1.8 ratio threshold.** With sequence truth available,
check the coordinate-only rule's sensitivity and specificity, and whether 1.8 is
the right cut. Currently 6 members deconvolved, 0 abstentions in the 12-cluster
batch.

---

## D. The clean re-run (needs the ORCID)

**Task D1 —** set `au_string` to the spelled-out name and add the ORCID prefix
in the Stockholm writer; confirm `stk lint` drops both the `au_format` ERROR and
the `orcid_missing` INFO.

**Task D2 —** re-run stages 0-2 end to end on the stratified batch, both merge
modes, both engines. Expect: 0 lint ERRORs, every alignment row verified against
the assembly, and `lint_triage.tsv` clean apart from anything genuinely new.

**Task D3 —** commit, and update `RESUME.md` headline numbers (they currently
describe the 12-cluster batch).

---

## E. Design decisions CG settled — implement, do not re-litigate

1. **Final classification comes from the curator.** `#=GF TP` as emitted is a
   *proposal*. The TIR-vs-ClassI structural cross-check is a **flag displayed
   beside the proposal**, never a gate, and the curator's classification
   overrides it in the deposited seed. (CG: *"the final label is changed
   according to curator decision."*)
2. **Structural / per-copy callers become a support channel, not clustering
   participants.** CG: *"split calls that are structural/per-copy vs homology
   (e.g. EDTA structural + LTRDeNovo) and include them when they overlap
   clusters of general purpose tools. They carry decent information about the TE
   being LTR and may drive good boundaries."*
   Implement `structural_ltr_support` per cluster: count overlapping fastltr
   members and ltrdenovo copies (ltrdenovo has one name per hit so it can never
   cluster — F5 — but every hit is classified LTR). Measured justification:
   fastltr joins 117/1,206 clusters, 96.6% of them LTR, and lifts the candidate
   rate 89% vs 56% — an effect that survives matching on which other tools are
   present. **Also record which member supplied the accepted boundary** once
   curator verdicts exist: CG's hypothesis is that structural callers have the
   boundaries right, and that is directly checkable.
3. **Linkage:** strict (mutual-best) is the default, lenient (components +
   density split at 0.5) runs alongside. B-lenient / mutual-top-2 was evaluated
   and rejected — it re-admits chaining.
4. **Merge mode:** gap_aware default, merge_always alongside for comparison.

---

## F. Answered questions — do not re-ask

- *Is the cluster-1103 element just a tandem repeat?* **No, measured.** Its
  6.87 Mb footprint overlaps fastan 1.16%, trf 0.56%, satellome 0.07% (G4
  rejects above 30%; all 12 built clusters fall between 0.006 and 0.089). It is
  an interspersed element that forms local tandem arrays — 2,892 same-strand
  runs, mostly of exactly two units. That combination is *why* its library
  consensus absorbed neighbouring units while the tandem finders stayed silent.
- *Can one genome settle the engine-disagreement question?* **Partly** — see §A.
  763 candidates give internal scale; cross-taxon generalisation and curator
  calibration do not follow from one assembly.

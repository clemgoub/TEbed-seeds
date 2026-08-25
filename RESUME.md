# Resume / handover — stages 0–1 working, step 3 next

Prototype status 2026-08-25. Design and measured rationale:
`PLAN_A_seedbuilder.md` rev 3 · algorithm deck `LHF_stage_algorithms.pptx` ·
**findings destined for the slides: `PIPELINE_FINDINGS.md`** (keep appending
there; the deck gets rebuilt from it).

---

## 0. Do these two things first

1. **`git init && git add -A && git commit`** — the repo is NOT under version
   control. The sandbox that scaffolded it could not run `git init` (broad
   host-grant restriction). Nothing here is committed.
2. **Spawn the parallel agent in §5** (consensus over-assembly test). It is
   independent of steps 1–4 and answers a question that could change how
   near-full-length is judged for a whole class of families.

## 1. How to run what exists

```bash
cd ~/Documents/TEbed-seeds
PYTHONPATH=. python -m lhfseeds.run_stage0 --linkage both   # ~80 s
PYTHONPATH=. python -m lhfseeds.run_stage1 --top 2          # ~30 s
```

Requires `pyyaml`, `pandas`, `numpy`, and the VGP_TEbed repo at the path in
`config/pipeline.yaml:vgp_tebed_repo`. Outputs go to `runs/<assembly>/`
(gitignored; `--work-dir` redirects).

## 2. Stage 0 — validated

Contract check → eligibility → equivalence graph → linkage → cluster gates.

GCA_951799975.1: 8,149 footprints, **6,625 edges**, **strict** 1,206 clusters /
max 14 / density 1.00 / **763 candidates**; **lenient** 1,081 / max 76 /
density 0.80 / 728 candidates.

Eligibility exclusions, both found by measurement and both load-bearing:
`windowmasker` (1 name for 4.36 M hits — one pseudo-family that links to
everything; caught by `mask_only`) and `ltrdenovo` (3,104 names / 3,104 hits =
per-copy ids; caught by `per_copy_ids`). Both still count for coverage/support
upstream — they just cannot cluster.

## 3. Stage 1 — working; three bugs found and fixed

1. **Consensus length must be MODAL, not max** (`stage1.consensus_length`).
   `rnd-1_family-286`: 31,387 hits imply 269 bp, 11 outliers imply 4,202 bp.
   Max made near-full-length read **0** vs a true ~14,600.
2. **Runaway merge chains capped.** Transitive merging walked tandem arrays
   into single copies at **329× the consensus** (32% of merge_always copies
   over 1.5×). Now split at the widest internal gap: max copy 88.5 kb →
   **14.4 kb**. See `PIPELINE_FINDINGS.md` F1 — this is a *slide-worthy*
   finding, not just a bugfix.
3. **The cap cannot run inside `merge_copies` for inherited lengths** (EDTA
   learns its length only after per-locus inheritance), so it is re-applied
   post-hoc in `run_stage1`.

Performance: the naive per-gap foreign-coverage scan was intractable (killed at
860 s having written nothing). Now pre-painted once per cluster: **7 s**
(`stage1.build_foreign_probe`).

### Per-locus inheritance, not cluster-wide

Cluster 62's mate consensus lengths are **bimodal** — 264/269 bp (rm2,
pantera) vs 764/765 bp (two REPET families). Any cluster-wide average is wrong
for both groups, so `build_clustermate_conslen` paints mate lengths per base
and each copy inherits from whichever mate covers its locus.
`cons_len_source ∈ bed16 | clustermate | none` (cluster 62: 114,019 / 45,280 /
2,113 rows). `is_full` stays **False** for inherited-length members rather than
guessed — near-full-length needs per-copy consensus coordinates, which EDTA
does not have. Flag, do not fake.

### Merge-mode comparison — first data on the "both ways" question

| cluster | mode | copies | near-full | multi-frag |
|---|---|---|---|---|
| 62 | merge_always | 83,692 | 29,136 | 881 |
| 62 | gap_aware | 84,369 | 28,845 | **214** |
| 540 | merge_always | 22,762 | 6,870 | 1,060 |
| 540 | gap_aware | 23,633 | 6,515 | **202** |

Gap-aware: more, shorter copies; multi-fragment copies down **4–5×**;
near-full-length essentially unchanged (−1%). Provisional default
**gap_aware** — but the decisive test is seed-level (depth, rebuilt consensus
length), i.e. step 3 below. **When you have it, update `PIPELINE_FINDINGS.md`
F3 and the stage-1.1 slide.**

---

## 4. Procedure to continue (step 3 onward)

Each step ends in a checkable state. Do them in order; 4.1 is the only one that
blocks on an external input.

**4.0 — commit** (§0.1), then create a branch: `git checkout -b stage2-seeds`.

**4.1 — get the assembly FASTA.** `config/pipeline.yaml:assembly_fasta` is
still `null`. Fetch GCA_951799975.1 (GenomeArk or NCBI, ~0.8 GB), set the path,
and verify sequence names match `data/GCA_951799975.1.chrom.sizes` in the hub
repo (the hub also ships a `chromAlias.txt` if they differ). *Check:*
`samtools faidx` succeeds and a random `copies.tsv` locus extracts non-N
sequence.

**4.2 — wire sampling into the runner.** `stage1.sample_copies` exists but is
not called. Add it after the merge step: cap 100 / floor 10 /
`full_length_frac` 0.60 from config. Record the **realized** composition (not
the target) per cluster — repeatome intactness varies by taxon, so the
achievable ratio does too. *Check:* every cluster emits ≤100 sampled copies and
`packet.json` carries realized counts.

**4.3 — extract sequence. Extract FRAGMENTS, never the merged span.**
`frag_starts`/`frag_ends` exist precisely so gap interiors are excluded — the
gap is often another element (62.9% measured). Add `flank_bp` (500) each side;
flanks are required for TSD detection and for Refiner to see where homology
stops. Strand: reverse-complement `-` copies. *Check:* extracted length equals
Σ fragment lengths + 2×flank, and a spot-check copy blasts back to its locus.

**4.4 — MAFFT consensus first** (fastest route to a lintable file; Refiner
comes later). Emit Stockholm: Smitten ids via `stage1.smitten` (already
converts BED 0-based half-open → 1-based fully closed), `.` gap character,
`#=GC RF` consensus line, required `#=GF DE/AU/TP/OC/SQ` (`au_string` from
config). `#=GF TP` needs `tp_map.tsv`, which does **not exist yet** — until it
does, emit seeds without `TP` and expect that specific lint failure.
*Check:* a `.stk` file exists for cluster 62 and parses.

**4.5 — `stk lint` and `stk lint --genome`** from the pinned `dfam-curator`
submodule (not yet added: `git submodule add https://github.com/Dfam-consortium/dfam-curator`).
**Soft gate** — a failing packet is still queued, flagged, with `lint.txt`
attached and its code counted in a batch `lint_triage.tsv`. *Check / project
milestone:* **one `stk lint`-clean seed built from track data.** That single
artifact de-risks the whole pipeline.

**4.6 — run 4.2–4.5 for both merge modes** on the same clusters and compare
seed depth and rebuilt consensus length. This is the decision point for
`merge_mode`; write the result into `PIPELINE_FINDINGS.md` F3.

**4.7 — Refiner vs MAFFT** on the same clusters; add the
`possible_overextension` guard (rebuilt consensus > 1.5× the **median** member
consensus length — median, not max, so one chimeric member cannot set the
reference).

**4.8 — `tp_map.tsv`** generation from the Dfam classification list, with
`scheme_version` and `source_url` columns. The Dfam scheme is being reconciled
with Repbase as the databases merge, so this must be cheap to regenerate and
every packet records the version it used. Unmapped canonical paths must **block**
seed emission rather than guess a `TP`.

Then: TE-Aid invocation (pooled + per-member tabs), queue UI, decisions with
classification and multi-reasons. `BRIEF_teaid_v2.md` is the TE-Aid v2 spec.

---

## 5. Parallel agent brief — consensus over-assembly test

**Spawn this as an independent agent; it does not touch stages 0–4.**

**Question.** Cluster 62 contains two REPET families whose consensi are 764 and
765 bp, while independent rm2 and pantera consensi for the same cluster are 269
and 264 bp. Coordinate-only evidence says the 765 bp entry is probably **~3
tandem units of a ~265 bp element collapsed into one library consensus**:

- REPET copies have median length 223 bp = **29% of their own consensus**
  (p75 = 266 bp = 35%);
- their match coordinates along the consensus cluster into **three blocks**
  (starts at 2 / 264 / 539; ends at 262 / 511 / 765);
- 765 / 3 ≈ 255 ≈ the rm2/pantera consensus length;
- rm2 hits form 2,892 same-strand tandem runs, median span 623 bp;
- the two REPET families partition loci almost perfectly (0.2% overlap of
  3,073 vs 3,280 copies) — the signature of redundant library entries.

A palindromic/MITE explanation was **tested and rejected**: adjacent rm2 hits
within 700 bp are 75% same-strand, only 25% opposite — below chance, whereas a
palindrome would show an excess.

**Tasks.**
1. Obtain the REPET consensus library FASTA (ask the user — he has said he will
   provide it at coding time) and the rm2/pantera library entries for
   `rnd-1_family-286` and `Unknown_572-fGobNig`.
2. **Self-align each 765 bp consensus** (blastn, `-word_size 7 -dust no`, or a
   dotplot). Three tandem blocks ⇒ off-diagonal same-strand repeats at ~±255
   and ~±510 bp offsets. Opposite-strand off-diagonals would instead support a
   palindrome and would overturn the coordinate-based reading.
3. **Align the 765 bp consensus against the 265 bp consensi.** Over-assembly
   predicts the short consensus matching the long one three times, at the block
   boundaries above.
4. Generalize into a **detector** that needs no library sequence, and test it
   across all 8,149 family footprints: flag a family when
   `median_copy_length ≈ consensus_length / k` for small integer k (2–5)
   **and** its match coordinates are block-structured rather than uniform.
   Report how many families are affected per tool — a per-program
   over-assembly rate is directly useful for library building.
5. Decide the pipeline consequence and write it up: when a member's consensus
   is an over-assembly, its own consensus length must **not** be used for the
   near-full-length judgement (it makes real copies look 71% truncated). Either
   inherit the unit length from cluster-mates (the machinery already exists:
   `build_clustermate_conslen`) or divide by the detected k. Add
   `cons_len_source = deconvolved` as a provenance value.

**Deliverables:** a short report with the alignment figures, the detector
implementation plus its per-tool rate table, and a patch (or a precise
recommendation) for how stage 1 should treat these families. Append the outcome
to `PIPELINE_FINDINGS.md` F4 — it is slated for a stage-1.5 slide, "when the
library consensus is the problem", which motivates the whole premise of
rebuilding consensi from copies.

---

## 6. Decisions already settled — do not re-litigate

- Linkage: **strict** default, lenient run alongside for comparison.
- Merge mode: **both**, decided at seed level (§4.6).
- Edge minimum **500 bp**, not 1 kb (measured: recovers MITE-scale links,
  94 edges, all involving families with median copy ≤ 800 bp).
- Order coherence ⅔, member-weighted, abstention-aware; **advisory per
  cluster**, never a filter.
- `stk lint` is a **soft** gate — failures must be visible, not silently
  dropped from the queue.
- Curator verdicts are `approved` / `rejected` plus `accepted_source`; there is
  no third "rebuild" verdict. `accepted_source ≠ pooled` triggers an automatic
  rebuild in the same batch.
- `au_string` is a placeholder pending group feedback.

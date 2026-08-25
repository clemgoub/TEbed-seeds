# Resume / handover — stages 0–2 working, lint-clean seeds built

Prototype status 2026-08-25 (session 2). Design and measured rationale:
`PLAN_A_seedbuilder.md` rev 3 · algorithm deck `LHF_stage_algorithms.pptx` ·
**findings destined for the slides: `PIPELINE_FINDINGS.md`** (keep appending
there; the deck gets rebuilt from it).

Branches: **`main`** is the prototype as it was handed over (one commit);
**`stage2-seeds`** carries this session's work and is where to continue.

---

## 0. The headline

**The project milestone is met: `stk lint`-clean Dfam seed alignments are being
built from track data, and they validate against the assembly.**

All seeds pass `stk lint --genome` with **ERROR=0, WARN=0**; the only remaining
diagnostic is the optional `orcid_missing` INFO. Independently:

- the rebuilt consensi reproduce RepeatModeler2's own library entries at
  **98.1%** (cluster 62: 265 bp vs 269 bp) and **95.2%** (cluster 540);
- **all 100** loci in RepeatModeler's own seed alignment for
  `rnd-1_family-286` are present in our 35,205-locus set for that cluster;
- MAFFT and Refiner, given the same 100 sampled copies, produce the **same
  265 bp consensus** (98.9% identical to each other).

---

## 1. How to run what exists

```bash
cd ~/Documents/TEbed-seeds
git checkout stage2-seeds
PYTHONPATH=. .venv/bin/python -m lhfseeds.run_stage0 --linkage both   # ~80 s
PYTHONPATH=. .venv/bin/python -m lhfseeds.run_stage1 --top 12         # ~13 min
PYTHONPATH=. .venv/bin/python -m lhfseeds.run_stage2 --top 12 \
      --modes gap_aware,merge_always --threads 4                      # ~6 min
PYTHONPATH=. .venv/bin/python -m pytest tests/ -q                     # 25 tests
```

`.venv/` is a Python 3.13 venv (`requirements.txt`; pandas is pinned `<3`).
Outputs go to `runs/<assembly>/` (gitignored; `--work-dir` redirects).

**Things now on disk that were missing before** (all gitignored, none in the repo):

| what | where | note |
|---|---|---|
| assembly FASTA | `data/GCA_951799975.1.fna` (+`.fai`) | was already in `~/Downloads`; 298/298 names **and** lengths match `chrom.sizes` |
| tool libraries | `libs/` | from GenomeArk `downstream_analyses/repeats/systematic_annotations/` |
| rm2's own seeds | `libs/GCA_951799975.1-families.stk` | 82 MB, RepeatModeler's Stockholm for this assembly — the natural comparison target |
| `stk` lint toolkit | `vendor/dfam-curator/` (submodule, Rust, `cargo build --release`) | binaries in `target/release/` |
| Refiner | `vendor/refiner/` | `tools/setup_refiner.sh` reproduces it |

`libs/` also has the **LTR_retriever** and **FastLTR** libraries — relevant to
the LTR thread in §5. The **REPET** library is still missing (you said you would
get it); `tools/overassembly_seqtest.sh` is written and waiting for it.

---

## 2. What was done this session (RESUME §4.0–4.8, all of it)

- **4.0** repo committed; work on branch `stage2-seeds`.
- **4.1** assembly FASTA sourced, indexed, verified.
- **4.2** sampling wired in — but see §3, pooled sampling needed a dedup rule
  the plan did not anticipate.
- **4.3** extraction: fragments (never the merged span), `flank_bp` both sides,
  reverse-complement on `-`. Two files per packet: `copies.fa` (element only,
  what gets aligned) and `copies.flanked.fa` (± 500 bp, for TSD/TE-Aid/Refiner).
- **4.4** MAFFT → consensus → Stockholm.
- **4.5** `stk lint` and `stk lint --genome`, soft gate, `lint.<engine>.txt`
  per packet plus a batch `lint_triage.tsv`.
- **4.6** both merge modes at seed level — **see F8, the answer is not what the
  plan expected**.
- **4.7** Refiner installed and wired as a second engine.
- **4.8** `config/tp_map.tsv` generated (`tools/build_tp_map.py`), TP resolved
  through the cluster's canonical path.

New code: `lhfseeds/fasta.py`, `stage2.py`, `run_stage2.py`, `stockholm.py`;
`tools/{build_tp_map,detect_overassembly}.py`, `tools/setup_refiner.sh`,
`tools/overassembly_seqtest.sh`; `tests/test_stage2.py`.

---

## 3. Things found that change the design

Full write-ups with numbers are in `PIPELINE_FINDINGS.md` F4 and F7–F11.
The ones that change how the pipeline works:

**Pooling must deduplicate by locus (F7).** A cluster's members are different
tools describing the same element, so they annotate the same loci over and
over: 90,313 copy rows in cluster 62 are 35,194 distinct loci, 2.57 annotations
each. Naive pooled sampling would have put the same sequence in the alignment
~2.6× over. Rule: one representative per locus, loci defined by **≥50%
reciprocal overlap** — adjacency is not enough, because this element is
tandemly arrayed and a touch-based rule chained whole arrays into single
"loci". That is F1's failure mode a third time, at a third level.

The payoff is that cross-tool support can then be counted, and it predicts copy
completeness: near-full-length fraction goes **0.4% → 25.7% → 63.6% → 73.4%**
for loci seen by 1, 2, 3, 4 tools. That is the first direct evidence that
multi-tool agreement carries information about copy *quality*, not just overlap.

**`stk lint --genome` earns its place (F9).** It caught two bugs that pass
tier-1 validation cleanly and would have shipped: MAFFT silently
reverse-complements rows and marks them `_R_` (**77 of 100 rows** in cluster
540) so the identifier's strand goes stale; and EDTA leaves the strand uncalled
(`.`, ~0.3% of its hits) which extraction and identifier-retagging were
resolving in opposite directions. Also: the placeholder `au_string`
`"C. Goubert, Wheeler Lab"` is a lint **ERROR** (`au_format` rejects
abbreviated first names) — it is now `"Clement Goubert"`, and an ORCID prefix
would clear the last INFO.

**Two hash-order bugs, and the second one moved a scientific value.** Python
set iteration depends on `PYTHONHASHSEED`, and this code accumulated through
sets in two places.

*Stage 0:* two runs over identical inputs produced the same clusters in a
different order, so `cluster_id` — a positional index — silently pointed at a
different family. Measured: cluster 62 became cluster 292 on a re-run, and id
62 came back holding an unrelated family, with no error raised. Every seed
packet keyed by id would have been invalidated. Clusters are now emitted in a
deterministic order (verified byte-identical across three hash seeds) and carry
a content-addressed `cluster_key`. **Use `cluster_key`, not `cluster_id`, to
refer to a cluster across runs** — the ids in older notes and in the first
draft of F4 are stale.

*Stage 1:* `build_clustermate_conslen` paints every mate's consensus length
into one shared array, so at a contested base whoever paints last wins — and
the caller passed a set. Measured across two runs with no other change: one
EDTA member's inherited consensus length moved from **381 bp to 1,755 bp**.
That is the number near-full-length is judged against. Members are now painted
in order of increasing evidence, so the best-supported mate wins a contested
base — deterministic, and the tie-break is a rule rather than an accident.

**Open design question this exposes:** when two cluster-mates genuinely
disagree about the element length at the same locus, "most hits wins" is a
defensible tie-break but it is a choice, not a derivation. Worth revisiting.

**An over-assembled member consensus wrecks the near-full-length judgement.**
Cluster 62's seed was admitting 495 and 642 bp copies of a ~265 bp element,
because REPET's consensus for them is 764 bp, so the 1.5× runaway cap was
1,146 bp. Two fixes: a copy is judged against the **cluster's modal** consensus
length in stage 2, and stage 1 now **deconvolves** a member whose own consensus
is ≥1.8× the median of its agreeing mates, recording
`cons_len_source = deconvolved`.

Effect on the 12-cluster batch: **6 members deconvolved, 0 abstentions**, and
near-full-length counts move by **24× to 494×** (F4 has the table). The largest
case is **rm2**, not REPET — `rnd-4_family-1870`, 1,992 bp against mates at
322 bp. Over-assembly is not a REPET quirk, and the first pass of this
investigation was wrong to conclude it was.

---

## 4. Procedure to continue

**4.1 — TE-Aid.** `~/Documents/TE-Aid` exists and `BRIEF_teaid_v2.md` is the
spec. Feed it `copies.flanked.fa` + `consensus.<engine>.fa` per packet, pooled
and per-member tabs.

**4.2 — the curator queue UI**, verdicts `approved`/`rejected` plus
`accepted_source`, with multi-reasons and classification.

**4.3 — scale past 12 clusters.** Stage 1 is the bottleneck at ~65 s/cluster
because `member_hits` and `build_clustermate_conslen` each re-read whole tool
BEDs per member. Reading each BED once into a per-family index would make the
763-candidate run tractable; nothing else needs to change.

**4.4 — per-member seeds.** Only pooled seeds are built. The design calls for
pooled + per-member tabs, and `accepted_source ≠ pooled` triggers a rebuild.

**4.5 — `possible_overextension` is computed but never acts.** It is recorded
in `packet.json` against the median member consensus length; decide whether it
should gate or only flag.

---

## 5. Open, and needing you

- **REPET library FASTA** — the last piece of the §5 sequence test.
  `tools/overassembly_seqtest.sh` runs the whole battery once it lands.
- **`au_string`** is `"Clement Goubert"` pending group feedback; add an ORCID
  to clear the last lint INFO.
- **`repeat:TE:ClassII:TIR:Academ` has no Dfam TP.** Dfam has `DNA/Academ-1`,
  `-2`, `-H` but no bare `DNA/Academ`, while every other TIR superfamily has a
  short form. 2 candidate clusters are affected; they emit without TP rather
  than assert a subfamily. Worth raising with Dfam.
- **The LTR thread you raised (F11).** fastltr joins only 117 of 1,206 clusters
  but **96.6% of them are LTR**, and its presence raises the candidate-gate
  pass rate from **56% to 89%** — an effect that survives matching on exactly
  which other tools are present. Structural LTR evidence looks like a *support
  channel* worth adding (`ltrdenovo`'s 3,104 per-copy calls are a second,
  currently unused one). Open question in F11: if "low-hanging fruit" means
  cheapest to curate confidently, fastltr-supported LTR clusters may be a
  better first batch than the highest-copy-number families.

---

## 6. Decisions settled — do not re-litigate

- Linkage: **strict** default, lenient run alongside for comparison.
- Merge mode: **gap_aware** — confirmed, but see F8: the seed-level test the
  plan expected to be decisive was *not* decisive (identical consensus length
  and depth). Gap-aware wins on the secondary criteria: more loci, fewer
  bridged rows, fewer over-length copies.
- Edge minimum **500 bp**, not 1 kb.
- Order coherence ⅔, member-weighted, abstention-aware; **advisory per
  cluster**, never a filter.
- `stk lint` is a **soft** gate — failures visible in `lint_triage.tsv`, never
  silently dropped from the queue.
- Curator verdicts are `approved` / `rejected` plus `accepted_source`.
- Seeds carry one alignment row **per fragment**, not per copy: each row is
  then a genuine contiguous interval, so gap interiors stay excluded *and*
  every Smitten identifier stays valid.
- An unmapped canonical path **blocks** `#=GF TP` rather than guessing it.

# Handoff — the LHF queue needs a size + protein-domain gate (2026-08-26)

**Prompted by CG:** *"the top 20 is contaminated by ultra small mites and other
likely non-autonomous things that are not acceptable for LHF — we need to
filter for size and consider prot domain (maybe expensive, but done once)."*

Correct, and measured. The current ranking (`candidates_strict.tsv` sorted by
`pooled_full_len`, `run_stage1.py:46`) is a **copy-number** ranking, so it
promotes exactly the wrong thing: abundant short non-autonomous elements.

## The contamination, quantified

Pipeline top-20 as ranked today: **median rebuilt consensus 469 bp; 19 of 20
under 1 kb; 12 under 500 bp**. That is MITE / SINE / fragment scale.

Across all 235 post-fix rebuilt models (213 clusters):

| tier | rule | n |
|---|---|---|
| A | core replication domain + class-appropriate size | **42** |
| B | core domain but short for its class | 4 |
| C | accessory domain only (likely truncated) | 7 |
| D | oversized for class, NO domain (suspect chimera/rDNA) | 15 |
| E | non-autonomous at expected size (SINE) | 2 |
| F | too short for class (MITE / fragment scale) | **165 (70%)** |

Note on the SINE rule: SINEs are exempt from the per-order minimum length (they
are non-autonomous by definition), which in the first pass let an *oversized*
SINE fall through to tier F ("too short for class") — nonsense for a 2926 bp
model. Fixed: a SINE above 700 bp now goes to tier D as oversized-for-class.
That moved one model (cluster 1161 `model1_2936bp`) from F to D and is the
difference between the counts above and the first version of this table.

## The filter (cheap — this is the key finding)

The protein evidence is **not** expensive. TE-Aid already ships 130 curated
TE-protein Pfam profiles (`dev-data/protein-cache/*/pfam.bhmm`, 7.3 MB).
Scanning **all 235 consensi (0.50 Mb) took 1.8 seconds** with
`bathsearch --cpu 8 -E 1e-5`. This can run inside stage 2 on every packet at
negligible cost — no need to treat it as a one-off.

Two-part gate as implemented (`lhf_size_domain_filter.tsv` has every model):

1. **Per-order minimum consensus length**: LTR/DIRS 4000, LINE 2500, RC 3000,
   DNA 1500, PLE 1200 bp; SINE exempt (non-autonomous by definition).
2. **At least one CORE replication domain at >=50% HMM coverage** — RT family
   (RVT_*, Pao_retrotransp, RT_RNaseH*, Exo_endo_phos*), transposase family
   (DDE_Tnp_*, Transposase_*, Dimer_Tnp_hAT, DBD_Tnp_*, MULE, Mutator,
   Helitron_like_N, Rol_Rep_N), or integrase (rve*, Integrase_H2C2, IN_DBD_C).

Requiring a *core* domain rather than *any* domain matters: a 461 bp PLE model
passed on `GIY-YIG` alone (endonuclease, no RT) — that is a fragment, and the
core-domain rule rejects it.

## Why both halves are needed

- Domain presence by size class: <300 bp **0%**, 300-600 **1.5%**,
  600-1200 **1.7%**, 1200-2500 **12%**, >2500 bp **84.5%**. Size alone is a
  decent proxy...
- ...but **5 short models carry a strong core domain** and would be wrongly
  discarded by a size cut alone — e.g. cluster 875, 1690 bp LINE/I-Jockey with
  `RVT_1` at 97% HMM coverage; cluster 848, 1294 bp with `RVT_1`.
- And **15 oversized models carry NO domain at all** — these are the suspects a
  size-only filter would wrongly *promote*. Largest is cluster 1091 `single` at
  13342 bp (order Unknown), then cluster 1146 `single` at 5464 bp (LTR/Pao,
  floor 4000) and cluster 890 `single` at 5328 bp. Cluster 1161 contributes
  **two** tier-D rows, and they are different models — do not conflate them:
  `model1_2936bp` (tp `SINE/5S`, 2926 bp; a 5S SINE should be ~120 bp, so this
  is almost certainly a chimera or rDNA array — it was the model shipped in the
  v4 demo queue) and `single` (2068 bp, order Unknown, no tp assigned).

## Open issues for the pipeline

1. **No DNA transposon reaches tier A.** Verified this is a property of the
   candidate set, not the profile library: the longest DNA model in the batch
   is 1577 bp and the longest RC/Helitron 1408 bp, while the profile set does
   contain DDE_Tnp_*, Dimer_Tnp_hAT, MULE, Mutator, Helitron_like_N. Either
   autonomous DNA elements are not surviving clustering/rebuilding, or this
   genome's DNA content really is non-autonomous. **Worth investigating** — a
   queue with zero DNA transposons is not a credible LHF deliverable.
2. **Add an alignment-depth floor.** Cluster 1169 passed tier A on a seed with
   **5 alignment rows at depth 3x** — domain-positive but not a trustworthy
   consensus. Suggest `aln_rows >= 10 AND median_depth >= 5` as a tier-A
   precondition. Flagged in the demo rather than dropped, for discussion.
3. **Ranking, not just filtering.** Within tier A, ranking by number of
   distinct core domains then full-length copies puts complete Gypsy elements
   (6 domains at 100% HMM coverage) at the top. That is a defensible LHF order;
   `pooled_full_len` is not.

Data: `lhf_size_domain_filter.tsv` (all 235 models, tier + domains),
`tierA_queue20.tsv` (the 20 shipped in the demo), `domains.tbl` (raw
bathsearch tblout).

# Handoff — rebuilt consensus regresses to the short length mode (CG review, 2026-08-26)

From CG's review of the curator demo (six top clusters from the 207-cluster
batch): *the rebuilt consensus is systematically no better than the best
member — often REPET's — and always lands near the shortest member (close to
RM2). Cluster 1183 has three clean full-length member consensi (REPET,
fastLTR reverse-complement, Pantera), yet the seed came out at 453 bp.*

This file is the measured diagnosis and the requested changes. Everything
below was measured on `runs/GCA_951799975.1` (gap_aware packets) in the
2026-08-26 session; re-derive before trusting numbers against a newer run.

## Diagnosis (measured on cluster 1183)

Member consensus lengths are BIMODAL — solo-LTR-sized vs full-element:

| member | cons_len |
|---|---|
| edta:TE_00000089 · rm2:ltr-1_family-35 · pantera:Unknown_466 | 423 / 423 / 414 |
| repet:G3000-Map20 · fastltr:CONS_4-7367 · pantera:Gypsy_9 · rm2:ltr-1_family-36 | 7400 / 7367 / 7359 / 6523 |

1. `run_stage2.py` takes `modal_cons = pool.cons_len.mode()` over COPY ROWS.
   Copy rows: 7400→89, 7367→45, **423→44**, 6523→13 over loci; but over the
   pre-dedup pool the 423 bp members dominate by row count → modal = 423.
   The mode is table-sensitive (pool mode 423 vs loci-table mode 7400) —
   fragile either way.
2. The F4 anti-tandem guard `max_copy_x_modal_consensus = 1.5` then drops
   every locus > 635 bp: **all 70 loci > 5 kb are excluded** (51 REPET,
   12 rm2:ltr-1_family-36, 6 fastltr, 1 edta) — `n_loci_over_modal = 83`.
3. Sampling therefore only ever sees the solo-LTR side (sampled max 618 bp);
   MAFFT rebuilds 453 bp, Refiner 426 bp. The "full-length" count in the
   packet (113) is full-length *relative to each member's own consensus*, so
   solo-LTR copies count as full and mask the loss.

The guard is doing what F4 intended (it was tuned on a 265 bp element whose
REPET entry was a ~3-unit tandem over-assembly); on a bimodal LTR cluster the
same rule deletes the full-length evidence instead of the artifact.

Batch scope: 8 / 207 clusters have max(member cons_len) > 2x modal —
LTR: 1183, 1201 · TIR: 5 clusters · LINE: 1. Small count, but these are
top-of-queue flagship clusters, so the visible damage is large.

## Requested changes (CG)

1. **Bimodality-aware modal length.** Before applying the x1.5 cap, test the
   member consensus lengths (or locus lengths) for bimodality (e.g. two
   modes > 2x apart). If bimodal:
   - prefer the LONG mode when it is supported by a structural LTR tool
     (fastltr / ltrdenovo member in the cluster is independent evidence the
     long form is a real full element, not a tandem multimer — see F11);
   - or build TWO consensi (solo-LTR + full element) and ship both in the
     packet — the LTR/INT pair is a legitimate two-model family in Dfam
     practice anyway.
   Compute the mode over DEDUPLICATED LOCI, not pool copy rows, so member
   row-count imbalance stops picking the mode (1183: loci mode is 7400).
2. **Do not let the F4 guard fire on loci that are near-full-length for a
   LONG member** — a locus spanning ~1.0x a 7.4 kb member consensus is not a
   tandem multimer of a 423 bp unit. The guard should compare against the
   length mode the locus itself belongs to.
3. **Refiner as the reported engine.** CG prefers Refiner. Flip
   `config/pipeline.yaml: consensus_engine` to `refiner` (both engines keep
   running; this changes which one is primary in packets.tsv and which seed
   the curator demo shows). F8b already found the engines agree at seed
   level, so this is presentation, not correctness.
4. **Re-run the affected clusters** (at least 1183, 1201 and the 5 TIR + 1
   LINE bimodal ones) and compare rebuilt_consensus_len against the best
   member; success = the rebuilt consensus reaches the full-element mode
   (~7.4 kb on 1183) or the packet ships both models explicitly.

## Pointers

- guard + modal: `lhfseeds/run_stage2.py` (~lines 158–180); F4 rationale in
  the comment there and in `PIPELINE_FINDINGS.md` F4.
- full-length definition (per-member, masks the loss): `lhfseeds/stage1.py`
  `near_full_length` + `run_stage1.py` (~line 254).
- structural-tool quality marker: `PIPELINE_FINDINGS.md` F11.
- EDTA library is on disk (`libs/GCA_951799975.1.fa.mod.EDTA.TElib.fa`) —
  note its entries for 1183 are headered `MITE/DTA`, which contradicts the
  cluster's LTR call; worth a look while in there.

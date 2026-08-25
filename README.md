# TEbed-seeds

Fast screen: candidate TE families from the VGP_TEbed multi-tool track →
consensus rebuilt from genomic copies → Dfam-compliant Stockholm seed →
TE-Aid QC → curator approve/reject with classification.

Design and measured rationale: `PLAN_A_seedbuilder.md` (rev 3); algorithm deck
`LHF_stage_algorithms.pptx`. **Measured findings live in
`PIPELINE_FINDINGS.md`** — the deck is rebuilt from that file. Current status
and what to do next: `RESUME.md`.

## Stages

| stage | module | does |
|---|---|---|
| 0 | `run_stage0` | eligibility → equivalence graph → strict/lenient linkage → cluster gates → candidates |
| 1 | `run_stage1` | gap-aware vs merge-always copy building, consensus-length inheritance |
| 2 | `run_stage2` | locus dedup → sampling → sequence extraction → MSA → consensus → Stockholm seed → `stk lint` |
| 3–4 | *not built* | TE-Aid QC, curator queue |

Stage 2 emits one seed packet per (cluster × merge mode), each carrying
`loci.tsv`, `sampled.tsv`, `copies.fa`, `copies.flanked.fa`,
`consensus.<engine>.fa`, `seed.<engine>.stk`, `lint.<engine>.txt` and a
`packet.json` recording realized (not target) composition and provenance.

## Quick start

```bash
python3.13 -m venv .venv && .venv/bin/pip install -r requirements.txt
git submodule update --init                 # vendor/dfam-curator
(cd vendor/dfam-curator && cargo build --release)
bash tools/setup_refiner.sh                 # optional, for the Refiner engine

PYTHONPATH=. .venv/bin/python -m lhfseeds.run_stage0 --linkage both
PYTHONPATH=. .venv/bin/python -m lhfseeds.run_stage1 --top 12
PYTHONPATH=. .venv/bin/python -m lhfseeds.run_stage2 --top 12 --threads 4
PYTHONPATH=. .venv/bin/python -m pytest tests/ -q
```

Set `assembly_fasta` in `config/pipeline.yaml`; stage 2 needs it. `runs/`,
`data/`, `libs/` and `vendor/refiner/` are gitignored.

## External requirements

`mafft`, `samtools`, `blastn` on PATH; a Rust toolchain to build `stk`; the
VGP_TEbed hub repo at `config/pipeline.yaml:vgp_tebed_repo`. Refiner
additionally needs RepeatMasker's Perl modules and a genuine `rmblastn` —
stock `blastn` cannot substitute; `tools/setup_refiner.sh` explains why and
installs both.

## Conventions worth knowing before editing

- BED is 0-based half-open; Dfam Smitten identifiers are **1-based fully
  closed** (`assembly:seq:start-end_strand`). Anything that trims or re-aligns
  a row must re-derive its coordinates — `stk lint --genome` is the check.
- Dfam's gap character is `.`; MAFFT writes `-`. Normalise on read.
- `cluster_id` is a positional index and only meaningful within one stage-0
  run. Use `cluster_key` (content-addressed) to refer to a cluster across runs.
- Consensus lengths are **modal**, never max; over-extension guards compare
  against the **median** member, never the max.

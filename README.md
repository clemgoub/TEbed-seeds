# TEbed-seeds

Fast screen: candidate TE families from the VGP_TEbed multi-tool track ->
consensus rebuilt from genomic copies -> Dfam-compliant Stockholm seed ->
TE-Aid QC -> curator approve/reject with classification.

Design and measured rationale: PLAN_A_seedbuilder.md (rev 3); algorithm deck
LHF_stage_algorithms.pptx.

Stages: 0 candidates (eligibility, equivalence graph, strict/lenient linkage,
gates) - 1 copies (gap-aware merge, extraction windows, sampling) - 2/3
consensus + seed + lint - 4 curator queue.

`runs/` is gitignored; `--work-dir` redirects it.

NOTE: repo not yet under git (sandbox limitation) -- run `git init && git add -A
&& git commit` on first Claude Code session.

Submodules (pinned later): dfam-curator (stk lint), TE-Aid v2.

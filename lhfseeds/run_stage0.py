"""Run stage 0 end to end.

    python -m lhfseeds.run_stage0 [--config config/pipeline.yaml] [--work-dir DIR]
                                  [--linkage strict|lenient|both]

Contract check happens FIRST (PLAN_A section 8 step 1): required inputs are
verified against the live VGP_TEbed repo before any computation.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd
import yaml

from . import stage0, gates


def check_contract(repo: Path, tools: list[str]) -> dict:
    """Verify the sync contract before running (blocking)."""
    problems = []
    for rel in ["config/tools.tsv", "docs/INPUT_FORMAT.md"]:
        if not (repo / rel).exists():
            problems.append(f"missing {rel}")
    sizes = list((repo / "data").glob("*.chrom.sizes"))
    if not sizes:
        problems.append("no chrom.sizes under data/")
    fams = repo / "report" / "data" / "families.parquet"
    if not fams.exists():
        problems.append("missing report/data/families.parquet (evidence)")
    present = [t for t in tools if (repo / "inputs" / f"{t}.bed").exists()]
    if len(present) < 2:
        problems.append(f"fewer than 2 tool BEDs present ({present})")
    if problems:
        sys.exit("[contract] FAILED:\n  " + "\n  ".join(problems))
    return dict(tools_present=present, chrom_sizes=str(sizes[0]))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    ap.add_argument("--work-dir", default=None)
    ap.add_argument("--linkage", default=None, choices=["strict", "lenient", "both"])
    args = ap.parse_args(argv)

    cfg = yaml.safe_load(open(args.config))
    repo = Path(cfg["vgp_tebed_repo"]).expanduser()
    work = Path(args.work_dir or cfg["work_dir"]).expanduser() / cfg["assembly"]
    work.mkdir(parents=True, exist_ok=True)
    linkage = args.linkage or cfg["linkage"]

    manifest = pd.read_csv(repo / "config" / "tools.tsv", sep="\t", comment="#")
    tools_manifest = manifest.tool_id.tolist()
    # family clustering is only meaningful for library/TE-scoped tools --
    # tandem finders and maskers have hub-scale footprints that chain everything
    FAMILY_SCOPES = {"general_homology", "structural_te", "ltr_only"}
    tool_scope = dict(zip(manifest.tool_id, manifest.scope))
    contract = check_contract(repo, tools_manifest)
    print(f"[contract] OK: {len(contract['tools_present'])} tool BEDs present", file=sys.stderr)

    # eligibility (four states, PLAN_A 3.1, + scope)
    status = stage0.assess_tools(repo, contract["tools_present"])
    elig = pd.DataFrame([vars(s) | dict(
        scope=tool_scope.get(s.tool, ""),
        clusterable=s.clusterable and tool_scope.get(s.tool, "") in FAMILY_SCOPES,
    ) for s in status.values()])
    elig.to_csv(work / "tool_eligibility.tsv", sep="\t", index=False)
    clusterable = elig[elig.clusterable].tool.tolist()
    print(f"[stage0] clusterable tools: {clusterable}", file=sys.stderr)

    # class map from the hub repo (source of truth for vocabulary)
    sys.path.insert(0, str(repo))
    import os
    cwd = os.getcwd()
    os.chdir(repo)
    from vgptrack.vocab import ClassMap
    cm = ClassMap.load()
    os.chdir(cwd)

    t0 = time.time()
    fam_hits = {t: stage0.load_families(repo, t, cm) for t in clusterable}
    fam_order, fam_path = {}, {}
    for t, df in fam_hits.items():
        for f, o in df.attrs["fam_order"].items():
            fam_order[f"{t}:{f}"] = o
        for f, pth in df.attrs["fam_path"].items():
            fam_path[f"{t}:{f}"] = pth
    sizes_file = repo / "data" / f"{cfg['assembly']}.chrom.sizes"
    sz = pd.read_csv(sizes_file, sep="\t", names=["chrom", "size"])
    chrom_sizes = dict(zip(sz.chrom, sz["size"]))

    edges, foot, rev = stage0.build_graph(
        fam_hits, chrom_sizes, cfg["edge_min_joint_bp"], cfg["edge_min_reciprocal"])
    pd.DataFrame(edges).to_parquet(work / "edges.parquet", index=False)

    # evidence join (from the hub's families table where available)
    fams = pd.read_parquet(repo / "report" / "data" / "families.parquet")
    fams["member"] = fams.tool + ":" + fams.family
    evidence = fams[["member"]].copy()
    evidence["n_full_len"] = fams.get("n_full_len")
    evidence["genomic_bp"] = fams.get("genomic_bp")
    evidence["cov_ge3_frac"] = fams.get("cov_ge3")
    evidence["div_median"] = fams.get("div_median")
    evidence["frac_tandemtool"] = fams.get("frac_tandemtool")

    results = {}
    modes = ["strict", "lenient"] if linkage == "both" else [linkage]
    for mode in modes:
        clusters = (stage0.link_strict(edges) if mode == "strict"
                    else stage0.link_lenient(edges, cfg["density_tau"]))
        diag = stage0.diagnostics(clusters, edges, fam_order)
        cand = gates.evaluate_clusters(clusters, evidence, fam_order, cfg,
                                       fam_path=fam_path)
        cand.to_csv(work / f"candidates_{mode}.tsv", sep="\t", index=False)
        results[mode] = dict(diag=diag, n_candidates=int(cand.candidate.sum()))
        print(f"[stage0:{mode}] clusters={diag['n_clusters']} max={diag['max_size']} "
              f"density={diag['median_density']:.2f} order-mixed={diag['frac_order_mixed']:.2f} "
              f"candidates={results[mode]['n_candidates']}", file=sys.stderr)

    json.dump(dict(config=cfg, contract=contract, linkage_run=modes,
                   results=results, elapsed_s=round(time.time() - t0, 1)),
              open(work / "stage0_summary.json", "w"), indent=1)
    print(f"[stage0] done in {time.time()-t0:.0f}s -> {work}", file=sys.stderr)


if __name__ == "__main__":
    main()

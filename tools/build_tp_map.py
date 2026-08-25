#!/usr/bin/env python3
"""Generate ``config/tp_map.tsv`` -- canonical class path -> Dfam ``#=GF TP``.

RESUME.md step 4.8.  A Dfam seed alignment must carry a ``#=GF TP`` value that
``stk lint`` recognises; ``stk lint``'s check is *exact set membership* against
the Dfam classification list (see ``vendor/dfam-curator/src/dfam/cache.rs``
``load_classification_tsv`` + ``lint.rs`` ``tier2_tp``).  A prefix of a valid
path is NOT valid.  This script therefore does two things:

1.  holds the curated correspondence between this project's canonical class
    vocabulary (``VGP_TEbed/config/class_map.tsv``, grammar
    ``repeat:{kind}:{class}:{order}:{superfamily}``) and Dfam class strings;
2.  **validates every proposed TP against the live Dfam list** and refuses to
    emit any value the list does not contain.  Nothing is guessed: a canonical
    path with no exact Dfam counterpart is left out of the file entirely, which
    (per the design decision in RESUME.md) blocks seed emission for that path.

Cheap to regenerate.  Dfam and Repbase are being reconciled, so both sides of
this map move.  Re-running the script re-reads the Dfam list, re-derives
``scheme_version`` from what it actually read, and re-reports coverage.

Sources
-------
Dfam list  https://www.dfam.org/releases/current/infrastructure/class_ns.tsv
           two columns: RepeatMasker-style short form TAB long semicolon path.
           Both forms are accepted by ``stk lint``; this file emits the short
           form in ``tp`` and carries the long form in ``tp_long``.
           Downloaded by ``vendor/dfam-curator/target/release/update-cache
           classifications`` into the stk cache (``$STK_CACHE_DIR``, else the
           platform cache dir + ``/stk``).
Project    ``VGP_TEbed/config/class_map.tsv`` (``vgptrack.vocab.ClassMap``).

Usage
-----
    .venv/bin/python tools/build_tp_map.py                    # default paths
    .venv/bin/python tools/build_tp_map.py --fetch            # refresh Dfam list first
    .venv/bin/python tools/build_tp_map.py --no-usage         # skip the BED scan
    .venv/bin/python tools/build_tp_map.py --coverage-out X.tsv

Exit status is 0 even when paths are unmapped -- unmapped is a legitimate,
reported state, not a failure.  ``--strict`` turns any unmapped *used* path
into exit 1 for CI.
"""
from __future__ import annotations

import argparse
import collections
import datetime as _dt
import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path

# --------------------------------------------------------------------------
# Revision of the correspondence table BELOW.  Bump whenever a row changes.
# The Dfam side is versioned separately, from the file that was actually read.
TPMAP_REV = "1"

DFAM_CLASS_URL = "https://www.dfam.org/releases/current/infrastructure/class_ns.tsv"
DFAM_CURRENT_INDEX = "https://www.dfam.org/releases/current/"

REPO = Path(__file__).resolve().parent.parent

# --------------------------------------------------------------------------
# THE CORRESPONDENCE TABLE
#
#   canonical_path -> (dfam_short_form, status, note)
#
# status:
#   exact        the Dfam node means the same thing as the canonical path.
#                Safe to emit automatically.
#   placeholder  Dfam has NO node at this depth, so the TP is Dfam's own
#                documented fallback (`Unknown` = Interspersed_Repeat;Unknown,
#                Dfam_Seeds.md "use a broad type such as
#                Interspersed_Repeat;Unknown as a temporary placeholder").
#                Information is LOST relative to the canonical path.  The seed
#                emitter must NOT use these silently -- treat a placeholder as
#                a block unless the caller explicitly opts in.
#
# A canonical path that appears in class_map.tsv but NOT here has no honest
# Dfam counterpart and is deliberately absent from the output: it blocks.
# --------------------------------------------------------------------------
TP_TABLE: dict[str, tuple[str, str, str]] = {
    # ---- Class I / LTR -----------------------------------------------------
    "repeat:TE:ClassI:LTR":        ("LTR", "exact",
        "Order-level; Long_Terminal_Repeat_Element."),
    "repeat:TE:ClassI:LTR:Gypsy":  ("LTR/Gypsy", "exact",
        "Ty3/Gypsy. Dfam is renaming Gypsy->Ty3; short form still LTR/Gypsy in this release."),
    "repeat:TE:ClassI:LTR:Copia":  ("LTR/Copia", "exact",
        "Dfam long form Ty1-Copia."),
    "repeat:TE:ClassI:LTR:BelPao": ("LTR/Pao", "exact",
        "Dfam short form is LTR/Pao; long form Bel-Pao."),
    "repeat:TE:ClassI:LTR:ERV":    ("LTR/ERV", "exact",
        "Retroviridae;Orthoretrovirinae. No ERV subgroup asserted."),
    "repeat:TE:ClassI:DIRS":       ("LTR/DIRS", "exact",
        "Wicker DIRS order = Dfam Tyrosine_Recombinase_Elements;DIRS."),
    # ---- Class I / non-LTR -------------------------------------------------
    "repeat:TE:ClassI:LINE":       ("LINE", "exact", "Order-level."),
    "repeat:TE:ClassI:LINE:L1":    ("LINE/L1", "exact", ""),
    "repeat:TE:ClassI:LINE:RTE":   ("LINE/RTE", "exact",
        "Dfam RTE-like clade (Group-II;Group-2;RTE-like)."),
    "repeat:TE:ClassI:LINE:R2":    ("LINE/R2", "exact", ""),
    "repeat:TE:ClassI:LINE:I":     ("LINE/I", "exact", ""),
    "repeat:TE:ClassI:LINE:Jockey": ("LINE/I-Jockey", "exact",
        "Dfam places Jockey inside the I-group: LINE/I-Jockey. Bare 'LINE/Jockey' is NOT in the Dfam list."),
    "repeat:TE:ClassI:SINE":       ("SINE", "exact", "Order-level; no promoter asserted."),
    "repeat:TE:ClassI:SINE:tRNA":  ("SINE/tRNA", "exact",
        "tRNA promoter, no core / LINE partner asserted."),
    "repeat:TE:ClassI:SINE:5S":    ("SINE/5S", "exact", "5S-RNA promoter."),
    "repeat:TE:ClassI:SINE:7SL":   ("SINE/7SL", "exact", "7SL-RNA promoter."),
    "repeat:TE:ClassI:PLE":        ("PLE", "exact", "Penelope-like elements."),
    # ---- Class II ----------------------------------------------------------
    "repeat:TE:ClassII:TIR":       ("DNA", "exact",
        "Wicker TIR order = Dfam Class_II_DNA_Transposition;Transposase, short form 'DNA'."),
    "repeat:TE:ClassII:TIR:hAT":          ("DNA/hAT", "exact", ""),
    "repeat:TE:ClassII:TIR:Tc1Mariner":   ("DNA/TcMar", "exact", "Dfam Tc1-Mariner."),
    "repeat:TE:ClassII:TIR:Mutator":      ("DNA/MULE", "exact", "Dfam Mutator-like."),
    "repeat:TE:ClassII:TIR:PIFHarbinger": ("DNA/PIF", "exact", "Dfam PIF-Harbinger."),
    "repeat:TE:ClassII:TIR:CACTA":        ("DNA/CMC", "exact",
        "CACTA superfamily = Dfam CACTA;CMC (CACTA-Mirage-Chapaev). Bare 'CACTA' is not a Dfam entry."),
    "repeat:TE:ClassII:TIR:PiggyBac":     ("DNA/PiggyBac", "exact", ""),
    "repeat:TE:ClassII:TIR:Kolobok":      ("DNA/Kolobok", "exact", ""),
    "repeat:TE:ClassII:TIR:Merlin":       ("DNA/Merlin", "exact", ""),
    "repeat:TE:ClassII:TIR:PElement":     ("DNA/P", "exact", "Dfam P_Element."),
    "repeat:TE:ClassII:TIR:Ginger":       ("DNA/Ginger", "exact", ""),
    "repeat:TE:ClassII:Helitron": ("RC", "exact",
        "Wicker Helitron order = Dfam Class_II_DNA_Transposition;Helicase, short form 'RC'. "
        "NOT RC/Helitron: that resolves to Helicase;Helitron-1 and would exclude Helentron/Helitron-2."),
    "repeat:TE:ClassII:Crypton":  ("DNA/Crypton", "exact",
        "Dfam Tyrosine_Recombinase;Crypton."),
    "repeat:TE:ClassII:Maverick": ("DNA/Maverick", "exact",
        "Dfam DNA_Polymerase;Maverick (Polinton)."),
    # ---- non-TE ------------------------------------------------------------
    "repeat:tandem:simple":    ("Simple_repeat", "exact", "Dfam Tandem_Repeat;Simple."),
    "repeat:tandem:satellite": ("Satellite", "exact",
        "Dfam Tandem_Repeat;Satellite; no centromeric/subtelomeric subtype asserted."),
    "repeat:multigene:tRNA":  ("tRNA", "exact", "Dfam Interspersed_Repeat;Pseudogene;RNA;tRNA."),
    "repeat:multigene:rRNA":  ("rRNA", "exact", ""),
    "repeat:multigene:snRNA": ("snRNA", "exact", ""),
    "repeat:multigene:scRNA": ("scRNA", "exact", ""),
    "repeat:artefact":        ("ARTEFACT", "exact", "Dfam 'Artifact'."),
    # ---- placeholders: Dfam has no node at this depth -----------------------
    "repeat": ("Unknown", "placeholder",
        "Abstention in this vocabulary (vgptrack UNINFORMATIVE_PATHS). Dfam Interspersed_Repeat;Unknown "
        "is the documented fallback but asserts 'interspersed', which 'repeat' does not."),
    "repeat:TE": ("Unknown", "placeholder",
        "Dfam has no 'transposable element, class unknown' node; the TE assertion is LOST."),
    "repeat:TE:ClassI": ("Unknown", "placeholder",
        "Dfam has no Class_I_Retrotransposition-only node; the Class I assertion is LOST."),
    "repeat:TE:ClassII": ("Unknown", "placeholder",
        "Dfam has no Class_II_DNA_Transposition-only node; the Class II assertion is LOST. "
        "Do NOT substitute 'DNA' -- that asserts the TIR order."),
}

# Canonical paths deliberately left unmapped, with the reason.  Reported, not
# emitted.  Keeping the reasons here makes regeneration self-documenting.
DELIBERATELY_UNMAPPED: dict[str, str] = {
    "repeat:TE:ClassII:TIR:Academ":
        "Dfam 4.0 has only DNA/Academ-1, -2, -H; no parent 'Academ' entry. "
        "Picking a numbered subgroup would be a guess.",
    "repeat:TE:ClassII:TIR:Sola":
        "Dfam 4.0 has only DNA/Sola-1, -2, -3; no parent 'Sola' entry.",
    "repeat:multigene:srpRNA":
        "'srpRNA' is a RepeatMasker class but is absent from Dfam class_ns.tsv in this release.",
    "repeat:tandem":
        "No bare 'Tandem_Repeat' entry; Simple_repeat and Satellite both assert more.",
    "repeat:multigene":
        "Ambiguous: covers RNA-gene and protein multigene families. Dfam 'RNA' "
        "(Interspersed_Repeat;Pseudogene;RNA) would assert RNA and 'interspersed'.",
}

BED16_COLS = [
    "chrom", "chromStart", "chromEnd", "name", "score", "strand", "SW_score",
    "perc_div", "perc_del", "perc_ins", "query_left", "repeat_class_family",
    "repeat_start", "repeat_end", "repeat_left", "hit_id",
]


# --------------------------------------------------------------------------
# Dfam classification list
# --------------------------------------------------------------------------
def stk_cache_dir() -> Path:
    """Mirror dfam-curator's cache_dir(): $STK_CACHE_DIR, else platform cache/stk."""
    env = os.environ.get("STK_CACHE_DIR")
    if env:
        return Path(env)
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Caches" / "stk"
    return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "stk"


def fetch_classification(dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["curl", "-fSL", "--retry", "2", DFAM_CLASS_URL, "-o", str(dest)],
                   check=True)


def load_classification(path: Path) -> tuple[set[str], dict[str, str]]:
    """Return (every accepted TP string, short_form -> long_form).

    Mirrors dfam-curator load_classification_tsv: BOTH columns are accepted, so
    membership in the returned set is exactly what `stk lint` will accept.
    """
    accepted: set[str] = set()
    long_of: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            t = line.strip()
            if not t or t.startswith("#"):
                continue
            if "\t" in t:
                short, long = (x.strip() for x in t.split("\t", 1))
                if short:
                    accepted.add(short)
                if long:
                    accepted.add(long)
                if short and long:
                    long_of[short] = long
            else:
                accepted.add(t)
    accepted.discard("rm_class")          # header row
    accepted.discard("dfam_lineage")
    long_of.pop("rm_class", None)
    return accepted, long_of


def release_from_previous(out: Path, sha: str) -> str:
    """Reuse the release recorded in a previous tp_map.tsv built from the same file.

    Lets an offline rebuild keep a real release number instead of 'unknown';
    only trusted when the Dfam file's sha256 is byte-identical.
    """
    if not out.exists():
        return "unknown"
    rel, seen_sha = "unknown", ""
    for line in out.read_text().splitlines():
        if not line.startswith("#"):
            break
        if line.startswith("#dfam_release:"):
            rel = line.split(":", 1)[1].strip()
        elif line.startswith("#dfam_class_sha256:"):
            seen_sha = line.split(":", 1)[1].strip()
    return rel if seen_sha == sha else "unknown"


def detect_release(offline: bool) -> str:
    """Read the RELEASE_DFAM_x_y marker file in the Dfam 'current' release dir."""
    if offline:
        return "unknown"
    try:
        out = subprocess.run(["curl", "-fsSL", "--max-time", "20", DFAM_CURRENT_INDEX],
                             capture_output=True, text=True, check=True).stdout
    except Exception as e:                                    # noqa: BLE001
        print(f"[tp_map] WARN cannot read {DFAM_CURRENT_INDEX}: {e}", file=sys.stderr)
        return "unknown"
    m = re.search(r"RELEASE_DFAM_(\d+)_(\d+)", out)
    return f"{m.group(1)}.{m.group(2)}" if m else "unknown"


# --------------------------------------------------------------------------
# Project vocabulary
# --------------------------------------------------------------------------
def load_class_map_vocab(hub: Path) -> tuple[set[str], str]:
    """Every canonical_path declared in the hub's class_map.tsv, + its version."""
    path = hub / "config" / "class_map.tsv"
    version, header, rows = "unknown", None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                body = line[1:].strip()
                if body.startswith("vocabulary_version:"):
                    version = body.split(":", 1)[1].strip()
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = [c.strip() for c in parts]
                continue
            rows.append(dict(zip(header, parts)))
    return {r["canonical_path"].strip() for r in rows if r.get("canonical_path", "").strip()}, version


def measure_usage(hub: Path, runs: Path, tools: list[str]) -> dict[str, dict]:
    """Canonical paths actually reached by this project's data.

    Family label -> canonical path exactly as lhfseeds.stage0.load_families does
    it (modal repeat_class_family per family name), then counted three ways:
    every family, families that landed in a cluster, and families in a cluster
    that passed the gates.
    """
    import pandas as pd
    sys.path.insert(0, str(hub))
    cwd = os.getcwd()
    os.chdir(hub)                       # ClassMap.load() uses a relative default
    try:
        from vgptrack.vocab import ClassMap
        cm = ClassMap.load()
    finally:
        os.chdir(cwd)

    fam_path: dict[str, str] = {}
    for tool in tools:
        bed = hub / "inputs" / f"{tool}.bed"
        if not bed.exists():
            print(f"[tp_map] WARN no BED for {tool}", file=sys.stderr)
            continue
        df = pd.read_csv(bed, sep="\t", names=BED16_COLS, skiprows=1,
                         usecols=["name", "repeat_class_family"], dtype=str)
        modal = df.groupby("name")["repeat_class_family"].agg(
            lambda x: x.mode().iat[0] if len(x.mode()) else "NA")
        for fam, lab in modal.items():
            fam_path[f"{tool}:{fam}"] = cm.lookup(str(lab), tool)["canonical_path"]

    clustered, candidate = set(), set()
    for f in sorted(runs.glob("candidates_*.tsv")):
        d = pd.read_csv(f, sep="\t")
        for _, r in d.iterrows():
            members = str(r["members"]).split(";")
            clustered.update(members)
            if bool(r["candidate"]):
                candidate.update(members)

    def tally(members) -> collections.Counter:
        return collections.Counter(fam_path[m] for m in members if m in fam_path)

    return {
        "all_families": collections.Counter(fam_path.values()),
        "clustered":    tally(clustered),
        "candidate":    tally(candidate),
        "n_families":   len(fam_path),
        "n_clustered":  len(clustered),
        "n_candidate":  len(candidate),
    }


# --------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--classification", default=None,
                    help="Dfam class_ns.tsv (default: stk cache classification.tsv)")
    ap.add_argument("--fetch", action="store_true",
                    help="download the Dfam list before building")
    ap.add_argument("--offline", action="store_true",
                    help="never touch the network (release detection returns 'unknown')")
    ap.add_argument("--hub", default=None, help="VGP_TEbed repo (default: from config/pipeline.yaml)")
    ap.add_argument("--runs", default=None, help="runs/<assembly> dir for the usage scan")
    ap.add_argument("--out", default=str(REPO / "config" / "tp_map.tsv"))
    ap.add_argument("--coverage-out", default=None, help="write the coverage table to this TSV")
    ap.add_argument("--no-usage", action="store_true", help="skip the tool-BED usage scan")
    ap.add_argument("--strict", action="store_true",
                    help="exit 1 if any canonical path used by the project is unmapped")
    args = ap.parse_args(argv)

    # ---- Dfam side -------------------------------------------------------
    cls_path = Path(args.classification) if args.classification \
        else stk_cache_dir() / "classification.tsv"
    if args.fetch or not cls_path.exists():
        if args.offline:
            sys.exit(f"[tp_map] {cls_path} missing and --offline given")
        print(f"[tp_map] fetching {DFAM_CLASS_URL}", file=sys.stderr)
        fetch_classification(cls_path)
    accepted, long_of = load_classification(cls_path)
    sha = hashlib.sha256(cls_path.read_bytes()).hexdigest()
    release = detect_release(args.offline)
    if release == "unknown":
        release = release_from_previous(Path(args.out), sha)
        if release != "unknown":
            print(f"[tp_map] reusing dfam_release {release} from the previous "
                  f"{args.out} (identical class_ns sha256)", file=sys.stderr)
    scheme_version = f"dfam-{release}+class_ns-{sha[:12]}+tpmap-r{TPMAP_REV}"
    print(f"[tp_map] {cls_path} -> {len(accepted)} accepted TP strings "
          f"(Dfam release {release}, sha256 {sha[:12]})", file=sys.stderr)

    # ---- validate every proposed TP against the live list ----------------
    rows, rejected = [], []
    for cpath, (tp, status, note) in sorted(TP_TABLE.items()):
        if tp not in accepted:
            rejected.append((cpath, tp))
            continue
        rows.append(dict(canonical_path=cpath, tp=tp, scheme_version=scheme_version,
                         source_url=DFAM_CLASS_URL, status=status,
                         tp_long=long_of.get(tp, ""), note=note))
    for cpath, tp in rejected:
        print(f"[tp_map] REJECTED {cpath} -> {tp!r}: not in the Dfam list; row dropped",
              file=sys.stderr)

    # ---- project side ----------------------------------------------------
    hub = Path(args.hub).expanduser() if args.hub else None
    if hub is None:
        import yaml
        cfg = yaml.safe_load(open(REPO / "config" / "pipeline.yaml"))
        hub = Path(cfg["vgp_tebed_repo"]).expanduser()
        runs = Path(args.runs).expanduser() if args.runs else \
            REPO / cfg["work_dir"] / cfg["assembly"]
    else:
        runs = Path(args.runs).expanduser() if args.runs else None

    vocab, vocab_version = load_class_map_vocab(hub)
    mapped = {r["canonical_path"] for r in rows}

    usage = None
    if not args.no_usage and runs and runs.exists():
        tools_tsv = hub / "config" / "tools.tsv"
        import pandas as pd
        man = pd.read_csv(tools_tsv, sep="\t", comment="#")
        fam_scopes = {"general_homology", "structural_te", "ltr_only"}
        elig = runs / "tool_eligibility.tsv"
        if elig.exists():
            e = pd.read_csv(elig, sep="\t")
            tools = e[e.clusterable].tool.tolist()
        else:
            tools = man[man.scope.isin(fam_scopes) & (man.ran.str.lower() == "yes")].tool_id.tolist()
        print(f"[tp_map] usage scan over {tools}", file=sys.stderr)
        usage = measure_usage(hub, runs, tools)

    # ---- write -----------------------------------------------------------
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    stamp = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    cols = ["canonical_path", "tp", "scheme_version", "source_url", "status", "tp_long", "note"]
    with open(out, "w") as fh:
        fh.write(f"""\
#tp_map_version: {scheme_version}
#generated: {stamp}
#generator: tools/build_tp_map.py (TPMAP_REV={TPMAP_REV})
#dfam_release: {release}
#dfam_class_url: {DFAM_CLASS_URL}
#dfam_class_sha256: {sha}
#dfam_class_file: {cls_path}
#dfam_accepted_strings: {len(accepted)}
#vocabulary_version: {vocab_version}
#vocabulary_source: {hub / 'config' / 'class_map.tsv'}
#
# canonical class path (vgptrack grammar repeat:kind:class:order:superfamily)
# -> Dfam '#=GF TP'.  `tp` is the RepeatMasker-style short form; `tp_long` is
# the equivalent semicolon path.  BOTH are accepted by `stk lint` (exact set
# membership -- a prefix of a valid path is NOT valid).
#
# status:
#   exact        the Dfam node means the same thing as the canonical path.
#   placeholder  Dfam has no node at this depth; `tp` is Dfam's documented
#                fallback and LOSES information.  The seed emitter must treat a
#                placeholder as a BLOCK unless the caller explicitly opts in.
#
# A canonical path absent from this file has no honest Dfam counterpart and
# MUST block seed emission (RESUME.md 4.8) -- never guess a TP.  Paths left out
# on purpose, with reasons, are listed at the bottom of this header.
#
# Every packet must record `tp_map_version` above.  Regenerate with:
#   .venv/bin/python tools/build_tp_map.py --fetch
#
# deliberately unmapped:
""")
        for cpath, why in sorted(DELIBERATELY_UNMAPPED.items()):
            fh.write(f"#   {cpath}\t{why}\n")
        fh.write("#\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")

    n_exact = sum(1 for r in rows if r["status"] == "exact")
    n_ph = len(rows) - n_exact
    print(f"[tp_map] wrote {out} -- {len(rows)} rows ({n_exact} exact, {n_ph} placeholder)",
          file=sys.stderr)

    # ---- coverage report -------------------------------------------------
    print(f"\n[coverage] class_map.tsv declares {len(vocab)} canonical paths "
          f"(vocabulary_version {vocab_version})", file=sys.stderr)
    missing_vocab = sorted(vocab - mapped)
    print(f"[coverage] mapped {len(vocab & mapped)}/{len(vocab)}; unmapped:", file=sys.stderr)
    for p in missing_vocab:
        why = DELIBERATELY_UNMAPPED.get(p, "NOT REVIEWED -- add to TP_TABLE or DELIBERATELY_UNMAPPED")
        print(f"    {p}\t{why}", file=sys.stderr)

    cov_rows = []
    unmapped_used = []
    if usage:
        print(f"\n[coverage] usage scan: {usage['n_families']} families, "
              f"{usage['n_clustered']} clustered members, "
              f"{usage['n_candidate']} candidate-cluster members", file=sys.stderr)
        for scope in ("all_families", "clustered", "candidate"):
            c = usage[scope]
            used = set(c)
            unm = sorted(used - mapped)
            tot = sum(c.values())
            covered = sum(v for k, v in c.items() if k in mapped)
            exact_only = sum(v for k, v in c.items()
                             if k in {r["canonical_path"] for r in rows if r["status"] == "exact"})
            print(f"[coverage] {scope}: {len(used - unm_set(unm))}/{len(used)} distinct paths mapped; "
                  f"{covered}/{tot} entries ({100*covered/max(tot,1):.1f}%); "
                  f"{exact_only}/{tot} via status=exact ({100*exact_only/max(tot,1):.1f}%)",
                  file=sys.stderr)
            for p in unm:
                print(f"    UNMAPPED {scope}: {p}  (n={c[p]})  "
                      f"{DELIBERATELY_UNMAPPED.get(p, 'NOT REVIEWED')}", file=sys.stderr)
            for p, n in sorted(c.items()):
                cov_rows.append(dict(scope=scope, canonical_path=p, n=n,
                                     status=next((r["status"] for r in rows
                                                  if r["canonical_path"] == p), "UNMAPPED")))
            if scope == "clustered":
                unmapped_used = unm

    if args.coverage_out and cov_rows:
        with open(args.coverage_out, "w") as fh:
            fh.write("scope\tcanonical_path\tn\tstatus\n")
            for r in cov_rows:
                fh.write(f"{r['scope']}\t{r['canonical_path']}\t{r['n']}\t{r['status']}\n")
        print(f"[tp_map] wrote {args.coverage_out}", file=sys.stderr)

    if args.strict and unmapped_used:
        return 1
    return 0


def unm_set(unm) -> set:
    return set(unm)


if __name__ == "__main__":
    sys.exit(main())

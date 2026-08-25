"""Stage 0 -- candidate identification.

Eligibility -> footprints -> equivalence graph -> linkage (strict|lenient)
-> cluster gates -> candidates.tsv + diagnostics.tsv.

Every design decision here was measured on GCA_951799975.1 before being coded
(PLAN_A_seedbuilder.md rev 3, section 0):
  - edge_min_joint_bp default 500: recovers MITE links (94 edges, 100%
    MITE-scale) that the earlier 1 kb rule starved.
  - strict (mutual-best) linkage default: chaining impossible by construction;
    connected components produced a 125-member, order-mixed cluster.
  - lenient = components + density-split at tau: retains more members
    (4,169 vs 3,599 measured) at density 0.81.
  - order-mixing is COMPARATIVE across linkages and advisory per cluster --
    it assumes tool labels, so it ranks strategies and informs curators;
    it never drops a copy.
"""
from __future__ import annotations

import collections
import hashlib
import itertools
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

BED16_COLS = [
    "chrom", "chromStart", "chromEnd", "name", "score", "strand", "SW_score",
    "perc_div", "perc_del", "perc_ins", "query_left", "repeat_class_family",
    "repeat_start", "repeat_end", "repeat_left", "hit_id",
]
NONFAMILY_PATH_PREFIXES = ("repeat:tandem", "repeat:multigene", "repeat:artefact")


# --------------------------------------------------------------------------- eligibility
@dataclass
class ToolStatus:
    tool: str
    ran: bool
    n_hits: int = 0
    n_names: int = 0
    frac_cons_coords: float = 0.0
    per_copy_ids: bool = False
    mask_only: bool = False

    @property
    def clusterable(self) -> bool:
        """Family clustering needs FAMILY identifiers: not one id per hit
        (ltrdenovo: 3,104 names / 3,104 hits) and not one id for everything
        (windowmasker: 1 name / 4.36M hits -- a single 'family' spanning
        hundreds of Mb links to every other family and chains the graph)."""
        return self.ran and not self.per_copy_ids and not self.mask_only

    @property
    def has_cons_coords(self) -> bool:
        return self.frac_cons_coords > 0.5


def assess_tools(repo: Path, tools: list[str]) -> dict[str, ToolStatus]:
    """Classify each tool: absent / ran+coords / ran-no-coords / per-copy ids.

    Absent tools leave every denominator (silence is never dissent).
    Per-copy-id tools (one name per hit) cannot cluster.
    """
    out: dict[str, ToolStatus] = {}
    for t in tools:
        bed = repo / "inputs" / f"{t}.bed"
        if not bed.exists():
            out[t] = ToolStatus(t, ran=False)
            continue
        n = 0
        names: set[str] = set()
        have_cc = 0
        with open(bed) as fh:
            next(fh)  # header
            for line in fh:
                fs = line.rstrip("\n").split("\t")
                n += 1
                names.add(fs[3])
                if fs[12] not in ("NA", "", "."):
                    have_cc += 1
        out[t] = ToolStatus(
            t, ran=True, n_hits=n, n_names=len(names),
            frac_cons_coords=have_cc / max(n, 1),
            per_copy_ids=(len(names) == n and n > 100),
            mask_only=(len(names) <= 3 and n > 1000),
        )
    return out


# --------------------------------------------------------------------------- footprints + graph
def load_families(repo: Path, tool: str, class_map) -> pd.DataFrame:
    """One tool's hits restricted to clusterable TE families."""
    df = pd.read_csv(
        repo / "inputs" / f"{tool}.bed", sep="\t", names=BED16_COLS, skiprows=1,
        usecols=["chrom", "chromStart", "chromEnd", "name", "repeat_class_family"],
    )
    modal = df.groupby("name")["repeat_class_family"].agg(
        lambda x: x.mode().iat[0] if len(x.mode()) else "NA")
    fam_path, fam_order = {}, {}
    for fam, lab in modal.items():
        p = class_map.lookup(str(lab), tool)["canonical_path"]
        parts = p.split(":")
        fam_path[fam] = p
        fam_order[fam] = parts[3] if len(parts) > 3 else (parts[1] if len(parts) > 1 else "repeat")
    keep = {f for f, p in fam_path.items()
            if not p.startswith(NONFAMILY_PATH_PREFIXES) and not f.startswith("(")}
    df = df[df.name.isin(keep)].copy()
    df.attrs["fam_order"] = {f: fam_order[f] for f in keep}
    df.attrs["fam_path"] = {f: fam_path[f] for f in keep}
    return df


def build_graph(fam_hits: dict[str, pd.DataFrame], chrom_sizes: dict[str, int],
                min_joint_bp: int, min_reciprocal: float, log=sys.stderr):
    """Per-chrom paint -> footprints and pairwise joint bp -> edges.

    Edge iff joint >= min_joint_bp AND joint >= min_reciprocal * smaller
    footprint. Weight channels: bp-reciprocal and copy co-annotation fraction.
    """
    tools = sorted(fam_hits)
    famid = {t: {f: i + 1 for i, f in enumerate(sorted(fam_hits[t].name.unique()))}
             for t in tools}
    rev = {t: {i: f for f, i in famid[t].items()} for t in tools}
    grp = {t: dict(tuple(fam_hits[t].groupby("chrom"))) for t in tools}

    foot: dict[tuple, int] = collections.defaultdict(int)
    joint = {p: collections.defaultdict(int)
             for p in itertools.combinations(tools, 2)}
    # copy co-annotation: copies of A whose midpoint lies in B's footprint
    copies_in = {p: collections.defaultdict(int)
                 for p in itertools.permutations(tools, 2)}
    n_copies: dict[tuple, int] = collections.defaultdict(int)

    for ch, L in chrom_sizes.items():
        arr = {}
        for t in tools:
            g = grp[t].get(ch)
            if g is None:
                arr[t] = None
                continue
            a = np.zeros(L, dtype=np.int32)
            fm = famid[t]
            for s, e, nme in zip(g.chromStart.to_numpy(), g.chromEnd.to_numpy(),
                                 g.name.to_numpy()):
                a[s:e] = fm[nme]
            arr[t] = a
            ids, cts = np.unique(a[a > 0], return_counts=True)
            for i, c in zip(ids, cts):
                foot[(t, int(i))] += int(c)
        for ta, tb in joint:
            A, B = arr[ta], arr[tb]
            if A is None or B is None:
                continue
            m = (A > 0) & (B > 0)
            if not m.any():
                continue
            combo = A[m].astype(np.int64) * (1 << 32) + B[m]
            ids, cts = np.unique(combo, return_counts=True)
            d = joint[(ta, tb)]
            for i, c in zip(ids, cts):
                d[(int(i >> 32), int(i & 0xFFFFFFFF))] += int(c)
        # copy midpoints
        for ta in tools:
            g = grp[ta].get(ch)
            if g is None:
                continue
            mid = ((g.chromStart.to_numpy() + g.chromEnd.to_numpy()) // 2)
            fa = g.name.map(famid[ta]).to_numpy()
            for f in np.unique(fa):
                n_copies[(ta, int(f))] += int((fa == f).sum())
            for tb in tools:
                if tb == ta or arr[tb] is None:
                    continue
                hit = arr[tb][mid]
                ok = hit > 0
                for f_a, f_b in zip(fa[ok], hit[ok]):
                    copies_in[(ta, tb)][(int(f_a), int(f_b))] += 1

    edges = []
    for (ta, tb), d in joint.items():
        for (ia, ib), c in d.items():
            fa, fb = foot[(ta, ia)], foot[(tb, ib)]
            if c >= min_joint_bp and c >= min_reciprocal * min(fa, fb):
                ca = copies_in[(ta, tb)].get((ia, ib), 0) / max(n_copies[(ta, ia)], 1)
                cb = copies_in[(tb, ta)].get((ib, ia), 0) / max(n_copies[(tb, ib)], 1)
                edges.append(dict(
                    a=f"{ta}:{rev[ta][ia]}", b=f"{tb}:{rev[tb][ib]}",
                    joint_bp=c, w_bp=c / min(fa, fb), w_copy=max(ca, cb)))
    print(f"[stage0] footprints={len(foot)} edges={len(edges)}", file=log)
    return edges, foot, rev


# --------------------------------------------------------------------------- linkage
def cluster_key(members) -> str:
    """Content-addressed cluster identifier.

    `cluster_id` is a positional index, so it is only meaningful within one
    stage-0 run: re-running with a different member set renumbers everything.
    The key is derived from the member set itself, so a seed packet can always
    be traced back to the cluster it was built from even across re-runs.
    """
    return hashlib.sha1(";".join(sorted(members)).encode()).hexdigest()[:10]


def _sorted_clusters(clusters):
    """Deterministic cluster order.

    Both linkages accumulate members through Python sets and dicts, whose
    iteration order depends on string hashing and therefore on PYTHONHASHSEED.
    Left alone, two runs of stage 0 over identical inputs produce identical
    CLUSTERS in a different ORDER, so `cluster_id` -- a positional index --
    silently points at a different family. Measured: cluster 62 became cluster
    292 on a re-run, and id 62 came back holding an unrelated family, which
    would have invalidated every seed packet keyed by it.
    """
    return sorted((sorted(c) for c in clusters), key=lambda c: (len(c), c))


def _components(pairs):
    par = {}
    def find(x):
        while par.setdefault(x, x) != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            par[ra] = rb
    comp = collections.defaultdict(list)
    for x in sorted(par):
        comp[find(x)].append(x)
    return _sorted_clusters(v for v in comp.values() if len(v) > 1)


def link_strict(edges):
    """Mutual-best per tool pair: chaining impossible by construction."""
    best = collections.defaultdict(dict)  # node -> other tool -> (w, partner)
    for e in edges:
        ta, tb = e["a"].split(":")[0], e["b"].split(":")[0]
        if best[e["a"]].get(tb, (-1, None))[0] < e["w_bp"]:
            best[e["a"]][tb] = (e["w_bp"], e["b"])
        if best[e["b"]].get(ta, (-1, None))[0] < e["w_bp"]:
            best[e["b"]][ta] = (e["w_bp"], e["a"])
    sel = {(n, m) for n, byt in best.items() for _, m in byt.values()}
    mutual = sorted((a, b) for a, b in sel if (b, a) in sel and a < b)
    return _components(mutual)


def link_lenient(edges, tau: float):
    """Connected components, then recursive density-split at the weakest cut."""
    E = collections.defaultdict(dict)
    for e in edges:
        E[e["a"]][e["b"]] = e["w_bp"]
        E[e["b"]][e["a"]] = e["w_bp"]

    def density(c):
        s = set(c)
        n = len(c)
        if n < 2:
            return 1.0
        ec = sum(1 for x in c for y in E[x] if y in s) / 2
        return ec / (n * (n - 1) / 2)

    def split(c, depth=0):
        if len(c) <= 2 or density(c) >= tau or depth > 12:
            return [c]
        s = set(c)
        intra = sorted((E[a][b], a, b) for a in c for b in E[a] if b in s and a < b)
        adj = {x: {y for y in E[x] if y in s} for x in c}
        for _, a, b in intra:
            adj[a].discard(b)
            adj[b].discard(a)
            seen, comps = set(), []
            for x in c:
                if x in seen:
                    continue
                stack, grpv = [x], []
                while stack:
                    y = stack.pop()
                    if y in seen:
                        continue
                    seen.add(y)
                    grpv.append(y)
                    stack.extend(sorted(adj[y] - seen))
                comps.append(sorted(grpv))
            if len(comps) > 1:
                out = []
                for g in comps:
                    if len(g) > 2:
                        out.extend(split(g, depth + 1))
                    elif len(g) > 1:
                        out.append(g)
                return out
        return [c]

    out = []
    for c in _components(sorted((e["a"], e["b"]) for e in edges)):
        out.extend(split(c))
    return _sorted_clusters(c for c in out if len(c) > 1)


def diagnostics(clusters, edges, fam_order):
    """The comparison table that justifies the linkage choice (per assembly)."""
    E = collections.defaultdict(set)
    for e in edges:
        E[e["a"]].add(e["b"])
        E[e["b"]].add(e["a"])
    sz = np.array([len(c) for c in clusters]) if clusters else np.array([0])
    dens, het = [], 0
    for c in clusters:
        s = set(c)
        n = len(c)
        ec = sum(1 for x in c for y in E[x] if y in s) / 2
        dens.append(ec / (n * (n - 1) / 2))
        orders = {fam_order.get(x) for x in c} - {None, "repeat"}
        if len(orders) > 1:
            het += 1
    return dict(
        n_clusters=len(clusters), members=int(sz.sum()), max_size=int(sz.max()),
        p99_size=int(np.quantile(sz, 0.99)) if len(sz) else 0,
        median_density=float(np.median(dens)) if dens else float("nan"),
        frac_order_mixed=het / max(len(clusters), 1),
    )

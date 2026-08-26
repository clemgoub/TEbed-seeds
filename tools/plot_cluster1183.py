"""Figure: why cluster 1183's consensus regressed to the short length mode.

Panel A  anatomy of one real full-element locus, showing that rm2 annotates the
         element in THREE pieces (LTR-INT-LTR) while fastltr, repet, pantera
         and edta each call it as ONE. That difference in annotation
         convention, not any disagreement about the element, is what makes the
         cluster's member consensus lengths bimodal.
Panel B  the consequence: the anti-tandem guard capped every locus at
         1.5 x 423 bp, deleting all 83 loci above 635 bp -- i.e. every
         full-length copy -- so the seed was rebuilt from solo LTRs alone.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
import pandas as pd, numpy as np

OUT = "docs/figures/d14_cluster1183_bimodal.png"
S = 39330767
LOC = "OX637597.1:39,330,768-39,338,142"

TRACKS = [   # label, [(rel_start, rel_end)], colour, note
    ("rm2:ltr-1_family-35", [(0, 423), (6952, 7375)], "#2b6cb0", "LTR model, 423 bp - both ends only"),
    ("rm2:ltr-1_family-36", [(423, 6952)], "#805ad5", "INT model, 6,523 bp - internal only"),
    ("fastltr:CONS_4", [(0, 7375)], "#2f855a", "whole element, 7,367 bp"),
    ("repet:G3000-Map20", [(-1, 7375)], "#2f855a", "whole element, 7,400 bp"),
    ("pantera:Gypsy_9", [(0, 7367)], "#2f855a", "whole element - a genuine member"),
    ("edta:TE_00000385", [(0, 7375)], "#718096", "whole element (no consensus coords)"),
]

fig = plt.figure(figsize=(13.5, 8.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1], hspace=0.42)

# ---------------------------------------------------------------- panel A
ax = fig.add_subplot(gs[0])
for i, (lab, segs, col, note) in enumerate(TRACKS):
    y = len(TRACKS) - i
    for (a, b) in segs:
        ax.add_patch(Rectangle((a, y - 0.28), b - a, 0.56, facecolor=col,
                               edgecolor="black", linewidth=0.6, alpha=0.9))
    ax.text(-260, y, lab, ha="right", va="center", fontsize=9.5, family="monospace")
    ax.text(7560, y, note, ha="left", va="center", fontsize=8.5, color="#4a5568")

# the element's own structure, drawn underneath
y0 = 0.15
ax.add_patch(Rectangle((0, y0 - 0.16), 423, 0.32, facecolor="#bee3f8",
                       edgecolor="#2b6cb0", linewidth=1.1))
ax.add_patch(Rectangle((6952, y0 - 0.16), 423, 0.32, facecolor="#bee3f8",
                       edgecolor="#2b6cb0", linewidth=1.1))
ax.add_patch(Rectangle((423, y0 - 0.10), 6529, 0.20, facecolor="#e9d8fd",
                       edgecolor="#805ad5", linewidth=1.1))
for x in (0, 423, 6952, 7375):
    ax.plot([x, x], [-0.35, 0.42], color="#a0aec0", lw=0.7, ls=":", zorder=0)
ax.text(211, y0 + 0.42, "LTR 423 bp", ha="center", fontsize=8.5, color="#2b6cb0")
ax.text(7163, y0 + 0.42, "LTR 423 bp", ha="center", fontsize=8.5, color="#2b6cb0")
ax.text(3688, y0 + 0.42, "internal region 6,529 bp   (one ORF, 1,898 aa: GAG zinc knuckle, RT YLDD)",
        ha="center", fontsize=8.5, color="#805ad5")
ax.text(3688, y0 - 0.62, "TG…CA termini · 4 bp TSD in 57 of 60 full-length loci · "
        "423 + 6,523 + 423 = 7,369", ha="center", fontsize=9, color="#2d3748")

ax.set_xlim(-2100, 10300); ax.set_ylim(-1.0, len(TRACKS) + 0.9)
ax.set_yticks([]); ax.set_xticks([0, 2000, 4000, 6000, 7375])
ax.set_xticklabels(["0", "2 kb", "4 kb", "6 kb", "7,375"])
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.set_title("A   One real locus, seven annotations — "
             f"{LOC}\n"
             "rm2 models the element in three pieces; the other four tools call it as one. "
             "That is the whole source of the bimodality.",
             fontsize=11.5, loc="left", pad=12)

# ---------------------------------------------------------------- panel B
ax2 = fig.add_subplot(gs[1])
loci = pd.read_csv("runs/GCA_951799975.1/seed_packets/cluster_01183/gap_aware/loci.tsv", sep="\t")
span = (loci.end - loci.start).to_numpy()
bins = np.logspace(np.log10(max(span.min(), 50)), np.log10(span.max() * 1.05), 60)
ax2.hist(span, bins=bins, color="#cbd5e0", edgecolor="#4a5568", linewidth=0.4)
ax2.set_xscale("log")
old_cap = 1.5 * 423
ax2.axvline(old_cap, color="#c53030", lw=2.2)
ax2.axvline(1.5 * 7363, color="#2f855a", lw=2.2, ls="--")
ax2.axvline(np.sqrt(418 * 7363), color="#4a5568", lw=1.2, ls=":")
n_lost = int((span > old_cap).sum())
ax2.text(old_cap * 1.06, ax2.get_ylim()[1] * 0.93,
         f"old cap 1.5 × 423 = {old_cap:.0f} bp\n→ {n_lost} loci deleted,\n   every full-length copy",
         color="#c53030", fontsize=9.5, va="top")
ax2.text(1.5 * 7363 * 0.97, ax2.get_ylim()[1] * 0.55, "new cap for the\nlong mode\n1.5 × 7,363",
         color="#2f855a", fontsize=9.5, va="top", ha="right")
ax2.text(np.sqrt(418 * 7363) * 1.05, ax2.get_ylim()[1] * 0.30,
         "mode boundary\n(geometric mean)", color="#4a5568", fontsize=8.5, va="top")
ax2.set_xlabel("locus span (bp, log scale)"); ax2.set_ylabel("loci")
for s in ("top", "right"): ax2.spines[s].set_visible(False)
ax2.set_title("B   Every locus was capped against the SHORT mode. "
              "Each locus is now capped against the mode it belongs to.",
              fontsize=11.5, loc="left", pad=10)

fig.savefig(OUT, dpi=170, bbox_inches="tight", facecolor="white")
print("wrote", OUT)
print(f"panel B: {len(span)} loci, {n_lost} above the old cap")

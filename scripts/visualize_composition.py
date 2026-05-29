#!/usr/bin/env python3
"""
visualize_composition.py -- tutorial input and result figures.

Produces:
  figures/references_by_class.png    horizontal bar chart of reference genomes by GTDB class
  figures/query_ground_truth.png     horizontal bar by source organism (novel highlighted)
  figures/dist_vs_truth.png          side-by-side: krepp dist --summarize vs ground truth

Usage:
  python3 scripts/visualize_composition.py \
      --metadata   data/metadata.tsv \
      --truth      data/query.truth.tsv \
      --dist-summ  precomputed/base/dist_summarize.tsv \
      --out-dir    figures/
"""

import argparse
import collections
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ---------------------------------------------------------------------------
# Color palettes
# ---------------------------------------------------------------------------

PHYLUM_COLORS = {
    "Cyanobacteriota":  "#2a9d8f",
    "Pseudomonadota":   "#e9c46a",
    "Bacteroidota":     "#f4a261",
    "Thermoproteota":   "#9b5de5",   # archaea
    "Bacillota":        "#adb5bd",   # outgroup
    "other":            "#ced4da",
}

NOVEL_COLOR = "#e63946"
TIP_DEFAULT = "#457b9d"


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

def read_metadata(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def read_truth(path):
    """Return {genome_id: read_count} and {genome_id: role} from truth TSV.
    Uses the 'count' column (total reads for that genome) to avoid double-counting.
    """
    counts = {}
    role_of = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            gid  = parts[idx["source_genome"]]
            role = parts[idx["role"]]
            n    = int(parts[idx["count"]])
            counts[gid] = n        # same value for all rows of this gid
            role_of[gid] = role
    return counts, role_of


def read_dist_summarize(path):
    """Return {reference_name: sequence_abundance} from krepp dist --summarize output."""
    abund = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("REFERENCE_NAME"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                abund[parts[0]] = float(parts[2])
    return abund


# ---------------------------------------------------------------------------
# Figure 1: reference genomes by class
# ---------------------------------------------------------------------------

def plot_class_distribution(rows, meta_by_id, out_path):
    classes = collections.Counter()
    for r in rows:
        if r.get("role") == "holdout":
            continue
        c = r.get("class") or "Unclassified"
        c = c.strip() or "Unclassified"
        classes[c] += 1

    items = sorted(classes.items(), key=lambda kv: (-kv[1], kv[0]))
    labels = [k for k, _ in items]
    values = [v for _, v in items]

    # color by phylum lookup
    phylum_of_class = {}
    for r in rows:
        phylum_of_class[r.get("class", "")] = r.get("phylum", "other")

    colors = []
    for lab in labels:
        ph = phylum_of_class.get(lab, "other")
        colors.append(PHYLUM_COLORS.get(ph, PHYLUM_COLORS["other"]))

    fig, ax = plt.subplots(figsize=(8, max(3.5, 0.55 * len(labels))))
    ax.barh(labels, values, color=colors, edgecolor="white", linewidth=0.6)
    ax.invert_yaxis()
    ax.set_xlabel("Number of reference genomes")
    ax.set_title(
        "Reference set by taxonomic class  (n={} genomes; holdout excluded)".format(sum(values)),
        fontsize=10,
    )
    for i, v in enumerate(values):
        ax.text(v + 0.15, i, str(v), va="center", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print("Wrote", out_path)


# ---------------------------------------------------------------------------
# Figure 2: query ground truth
# ---------------------------------------------------------------------------

def plot_query_ground_truth(counts, role_of, metadata, out_path):
    if not counts:
        print("No truth data; skipping query_ground_truth.png")
        return

    species_of = {r["genome"]: r.get("species", r["genome"]) for r in metadata}
    phylum_of  = {r["genome"]: r.get("phylum", "other")      for r in metadata}

    items  = sorted(counts.items(), key=lambda kv: -kv[1])
    labels = []
    sizes  = []
    colors = []
    for gid, n in items:
        sp = species_of.get(gid, gid)
        labels.append("{} ({})".format(sp, gid))
        sizes.append(n)
        if role_of.get(gid) == "novel":
            colors.append(NOVEL_COLOR)
        else:
            ph = phylum_of.get(gid, "other")
            colors.append(PHYLUM_COLORS.get(ph, TIP_DEFAULT))

    total = sum(sizes)
    fig, ax = plt.subplots(figsize=(10, max(4, 0.55 * len(labels))))
    ax.barh(range(len(labels)), sizes, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel("Read count")
    ax.set_title(
        "Simulated community ({} reads)  --  red = novel (held-out) organism".format(total),
        fontsize=10,
    )
    for i, n in enumerate(sizes):
        ax.text(n + total * 0.003, i, "{:.1f}%".format(100 * n / total),
                va="center", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend for phyla
    seen_ph = {}
    for gid, _ in items:
        if role_of.get(gid) == "novel":
            seen_ph["novel (held-out)"] = NOVEL_COLOR
        else:
            ph = phylum_of.get(gid, "other")
            seen_ph[ph] = PHYLUM_COLORS.get(ph, TIP_DEFAULT)
    patches = [mpatches.Patch(color=c, label=l) for l, c in seen_ph.items()]
    ax.legend(handles=patches, fontsize=8, loc="lower right", framealpha=0.7)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print("Wrote", out_path)


# ---------------------------------------------------------------------------
# Figure 3: krepp dist --summarize vs ground truth
# ---------------------------------------------------------------------------

def plot_dist_vs_truth(dist_abund, truth_counts, role_of, metadata, out_path):
    if not dist_abund or not truth_counts:
        print("Missing data for dist_vs_truth; skipping.")
        return

    species_of = {r["genome"]: r.get("species", r["genome"]) for r in metadata}
    total_truth = sum(truth_counts.values())

    # Only include genomes that appear in krepp dist output with > 0.1% abundance
    # or are in the ground truth with > 0.1%
    all_gids = set()
    for gid, n in truth_counts.items():
        if n / total_truth >= 0.001:
            all_gids.add(gid)
    for gid, ab in dist_abund.items():
        if ab >= 0.001:
            all_gids.add(gid)

    # Sort by ground truth abundance (NaN = 0 for genomes only in dist output)
    def truth_ab(gid):
        return truth_counts.get(gid, 0) / total_truth

    items = sorted(all_gids, key=truth_ab, reverse=True)
    n = len(items)
    labels  = [species_of.get(gid, gid) + "\n(" + gid + ")" for gid in items]
    truth_v = [truth_ab(gid) * 100 for gid in items]
    dist_v  = [dist_abund.get(gid, 0.0) * 100 for gid in items]

    x = range(n)
    bar_w = 0.38
    fig, ax = plt.subplots(figsize=(max(10, n * 0.72), 5))
    bars1 = ax.bar([i - bar_w / 2 for i in x], truth_v, bar_w,
                   label="Ground truth", color="#457b9d", edgecolor="white", linewidth=0.5)
    bars2 = ax.bar([i + bar_w / 2 for i in x], dist_v, bar_w,
                   label="krepp dist --summarize", color="#e9c46a",
                   edgecolor="white", linewidth=0.5)

    # Highlight novel bars
    for i, gid in enumerate(items):
        if role_of.get(gid) == "novel":
            ax.bar(i - bar_w / 2, truth_v[i], bar_w,
                   color=NOVEL_COLOR, edgecolor="white", linewidth=0.5)
            ax.annotate("novel\n(held-out)", xy=(i, max(truth_v[i], dist_v[i])),
                        ha="center", va="bottom", fontsize=7.5,
                        color=NOVEL_COLOR, fontweight="bold")

    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=7, rotation=30, ha="right")
    ax.set_ylabel("Relative abundance (%)")
    ax.set_title(
        "krepp dist --summarize vs. ground truth\n"
        "(novel PROC_AS9601 reads inflate PROC_MIT9301 estimated abundance)",
        fontsize=10,
    )
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print("Wrote", out_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata",  required=True)
    ap.add_argument("--truth",     required=True)
    ap.add_argument("--dist-summ", default=None,
                    help="dist_summarize.tsv from krepp dist --summarize (e.g. precomputed/base)")
    ap.add_argument("--out-dir",   required=True)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    rows       = read_metadata(args.metadata)
    meta_by_id = {r["genome"]: r for r in rows}

    plot_class_distribution(rows, meta_by_id, out / "references_by_class.png")

    counts, role_of = read_truth(args.truth)
    plot_query_ground_truth(counts, role_of, rows, out / "query_ground_truth.png")

    if args.dist_summ and Path(args.dist_summ).exists():
        dist_abund = read_dist_summarize(args.dist_summ)
        plot_dist_vs_truth(dist_abund, counts, role_of, rows, out / "dist_vs_truth.png")
    else:
        print("--dist-summ not provided; skipping dist_vs_truth.png")


if __name__ == "__main__":
    main()

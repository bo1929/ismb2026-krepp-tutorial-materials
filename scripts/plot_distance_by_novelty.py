#!/usr/bin/env python3
"""Plot krepp distance distributions stratified by novelty level."""

from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
QUERY_INFO = ROOT / "data" / "query_info.tsv"
QUERY_GENOMES = ROOT / "data" / "query_genomes"
RESULTS = ROOT.parent / "ismb2026-krepp-tutorial-materials" / "results" / "distances_filtered.tsv"
OUT_FIGURE = ROOT / "figures" / "distance_by_novelty.png"

# ── layout ──────────────────────────────────────────────────────────────────
FIG_DPI = 180
FIGSIZE = (8.0, 5.2)
BASE_FS = 11.5
SMALL_FS = 9.5
TITLE_FS = 13.0

NOVELTY_ORDER = ["species", "genus", "family", "kingdom"]
NOVELTY_LABELS = {n: n.capitalize() for n in NOVELTY_ORDER}
NOVELTY_COLORS = {
    "species": "#2f7d6c",
    "genus":   "#557fa5",
    "family":  "#a86b2f",
    "kingdom": "#aa4a46",
}
NOVELTY_PALETTE = {
    "species": "#43a88d",
    "genus":   "#6b9cc2",
    "family":  "#cc8a4a",
    "kingdom": "#c4625e",
}

XLABEL = "Novelty level of the query genome"
YLABEL = "Estimated distance (krepp)"
TITLE = "Per-read distance distributions by novelty level"


# ═══════════════════════════════════════════════════════════════════════════════
# Data loading
# ═══════════════════════════════════════════════════════════════════════════════

def load_query_info():
    """Return {taxid: {organism, novelty_level, ...}}."""
    by_taxid = {}
    with open(QUERY_INFO, encoding="utf-8") as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            tid = (rec.get("taxid") or "").strip()
            if tid:
                by_taxid[tid] = rec
    return by_taxid


def build_accession_map(query_by_taxid):
    """Map every FASTA header accession to its organism and novelty level."""
    acc_to_org = {}
    acc_to_novelty = {}
    for fasta in sorted(QUERY_GENOMES.glob("*.fna.gz")):
        taxid = fasta.stem.replace(".fna", "")
        q = query_by_taxid.get(taxid, {})
        org = q.get("organism", taxid)
        novelty = q.get("novelty_level") or "other"
        with gzip.open(fasta, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    acc = line[1:].split()[0]
                    acc_to_org[acc] = org
                    acc_to_novelty[acc] = novelty
    return acc_to_org, acc_to_novelty


def load_distances_and_mapping_rate(results_path, acc_to_org, acc_to_novelty):
    """Parse the filtered TSV.

    Returns
    -------
    rows : list[dict]
        Non-NaN rows with {org, novelty, ref, dist}.
    mapping_rate : dict[str, float]
        Fraction of unique reads per novelty level that got ≥1 non-NaN hit.
    """
    rows = []
    # Per novelty level: set of read IDs with ≥1 non-NaN hit, and total unique reads
    mapped_ids = defaultdict(set)
    all_ids = defaultdict(set)

    with open(results_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3 or parts[0] == "SEQ_ID":
                continue
            seq_id, ref, dist_s = parts[0], parts[1], parts[2]
            acc = seq_id.rsplit("-", 1)[0] if "-" in seq_id else seq_id
            org = acc_to_org.get(acc)
            novelty = acc_to_novelty.get(acc, "other")
            if org is None or novelty == "other":
                continue

            all_ids[novelty].add(seq_id)

            if dist_s != "NaN":
                dist = float(dist_s)
                rows.append({"org": org, "novelty": novelty, "ref": ref, "dist": dist})
                mapped_ids[novelty].add(seq_id)

    mapping_rate = {}
    for nov in NOVELTY_ORDER:
        total = len(all_ids.get(nov, set()))
        mapped = len(mapped_ids.get(nov, set()))
        mapping_rate[nov] = mapped / total if total > 0 else 0.0

    return rows, mapping_rate


# ═══════════════════════════════════════════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════════════════════════════════════════

def plot_distance_by_novelty(novelty_dists, mapping_rate, out_path):
    """Violin + strip plot of distances grouped by novelty level."""
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": BASE_FS,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titleweight": "bold",
    })

    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=FIG_DPI)

    levels = [n for n in NOVELTY_ORDER if n in novelty_dists]
    data = [novelty_dists[n] for n in levels]
    positions = list(range(len(levels)))

    # Violin plots
    vp = ax.violinplot(
        data, positions=positions,
        showmeans=False, showmedians=True, showextrema=False,
        widths=0.72,
    )
    for i, body in enumerate(vp["bodies"]):
        body.set_facecolor(NOVELTY_PALETTE[levels[i]])
        body.set_alpha(0.55)
        body.set_edgecolor(NOVELTY_COLORS[levels[i]])
        body.set_linewidth(1.0)
    vp["cmedians"].set_color("#2c2d34")
    vp["cmedians"].set_linewidth(2.2)

    # Overlay jittered strip points (subsample for performance)
    rng = np.random.default_rng(42)
    for i, d in enumerate(data):
        n_show = min(len(d), 600)
        sample = rng.choice(d, size=n_show, replace=False)
        jitter = rng.uniform(-0.22, 0.22, size=n_show)
        ax.scatter(
            [positions[i]] * n_show + jitter, sample,
            s=4.5, c=NOVELTY_COLORS[levels[i]], alpha=0.28,
            edgecolors="none", zorder=3,
        )

    # Stats annotation: mean distance + mapping rate
    for i, d in enumerate(data):
        arr = np.array(d)
        mean_v = np.mean(arr)
        mr = mapping_rate.get(levels[i], 0.0)
        y_pos = np.max(arr) + 0.018
        ax.text(
            positions[i], y_pos,
            f"μ = {mean_v:.3f}\nmapping rate = {mr:.0%}",
            ha="center", va="bottom", fontsize=SMALL_FS + 0.5,
            color="#5c5e69", linespacing=1.3,
        )

    # Axes
    ax.set_xticks(positions)
    ax.set_xticklabels([NOVELTY_LABELS.get(n, n) for n in levels], fontsize=BASE_FS)
    ax.set_ylabel(YLABEL, fontsize=BASE_FS + 1)
    ax.set_xlabel(XLABEL, fontsize=BASE_FS + 1, labelpad=8)
    ax.set_title(TITLE, loc="left", fontsize=TITLE_FS, fontweight="bold", pad=10)
    ax.set_ylim(-0.015, 0.255)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(0.05))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.tick_params(axis="y", labelsize=SMALL_FS)

    # Light horizontal reference lines
    for y in [0.05, 0.10, 0.15, 0.20]:
        ax.axhline(y=y, color="#d5d7de", linewidth=0.5, zorder=0)

    fig.tight_layout(pad=1.2)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight", pad_inches=0.06)
    fig.savefig(str(out_path).replace(".png", ".pdf"), dpi=FIG_DPI, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print(f"Wrote {out_path}")


# ═══════════════════════════════════════════════════════════════════════════════
def main():
    if not RESULTS.exists():
        raise FileNotFoundError(
            f"Results file not found: {RESULTS}\n"
            "Run krepp dist --filter first to produce distances_filtered.tsv"
        )

    query_by_taxid = load_query_info()
    acc_to_org, acc_to_novelty = build_accession_map(query_by_taxid)
    rows, mapping_rate = load_distances_and_mapping_rate(
        RESULTS, acc_to_org, acc_to_novelty,
    )

    # Group distances by novelty level
    novelty_dists = defaultdict(list)
    for r in rows:
        novelty_dists[r["novelty"]].append(r["dist"])

    # Print summary
    print("Distances by novelty level:")
    for nov in NOVELTY_ORDER:
        if nov in novelty_dists:
            arr = np.array(novelty_dists[nov])
            mr = mapping_rate.get(nov, 0.0)
            print(f"  {nov}: n={len(arr):,}  mean={arr.mean():.4f}  "
                  f"median={np.median(arr):.4f}  "
                  f"IQR=[{np.percentile(arr,25):.4f}, {np.percentile(arr,75):.4f}]  "
                  f"mapped={mr:.1%}")

    plot_distance_by_novelty(novelty_dists, mapping_rate, OUT_FIGURE)


if __name__ == "__main__":
    main()

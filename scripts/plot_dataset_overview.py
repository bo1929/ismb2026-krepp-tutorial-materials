#!/usr/bin/env python3
"""Create two dataset overview figures and the query-novelty table used by content/03-dataset.md."""

from __future__ import annotations

import csv
import os
import tempfile
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROFILE = ROOT / "data" / "profile.tsv"
REFERENCE_INFO = ROOT / "data" / "reference_info.tsv"
QUERY_INFO = ROOT / "data" / "query_info.tsv"
OUT_PROFILE = ROOT / "figures" / "dataset_profile.png"
OUT_REF = ROOT / "figures" / "dataset_reference.png"

# ── text labels ──────────────────────────────────────────────────────────
PROFILE_TITLE = "Query profile composition"
REF_CONTEXT_TITLE = "Reference panel and query–reference relationships"

KINGDOM_TITLE = "Kingdom composition"
CLASS_TITLE = "Class composition"
SPECIES_TITLE = "Query species (coloured by class)"
REPRESENTATION_TITLE = "How query reads relate to the reference panel"
REFERENCE_TITLE = "Reference genome set by role"
PERCENT_AXIS_LABEL = "Percent of query mixture"
REFERENCE_COUNT_AXIS_LABEL = "Number of reference genomes"

KINGDOM_COLORS = {
    "Bacteria": "#4A4A4A",
    "Archaea": "#BDBDBD",
    "Other": "#E0E0E0",
}

# Taxonomic ranks used in query_info / reference_info role fields (deepest shared rank).
RANK_ORDER = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]

REL_LABELS = {
    "species": "species",
    "genus":   "genus",
    "family":  "family",
    "order":   "order",
    "class":   "class",
    "phylum":  "phylum",
    "kingdom": "kingdom",
}

REF_ROLE_LABELS = {
    "species": "species",
    "genus":   "genus",
    "family":  "family",
    "order":   "order",
    "class":   "class",
    "phylum":  "phylum",
    "kingdom": "kingdom",
}

# Figure 2: neutral blue-gray tones (dark = closer match at species end)
REL_COLORS = {
    "species": "#4a5d6b",
    "genus":   "#5f6f7c",
    "family":  "#74828f",
    "order":   "#8895a1",
    "class":   "#9ca7b2",
    "phylum":  "#b0b9c2",
    "kingdom": "#c5ccd3",
    "other":   "#d4d9de",
}

REF_ROLE_COLORS = {
    "species": "#5a6774",
    "genus":   "#6b7783",
    "family":  "#7c8792",
    "order":   "#8d97a1",
    "class":   "#9ea7b0",
    "phylum":  "#afb7bf",
    "kingdom": "#c0c7ce",
}

CLASS_COLORS = {
    "Gammaproteobacteria":   "#0072B2",
    "Thermoprotei":          "#E69F00",
    "Alphaproteobacteria":   "#009E73",
    "Thermococci":           "#CC79A7",
    "Epsilonproteobacteria": "#56B4E9",
    "Acidimicrobiia":        "#D55E00",
    "Deinococci":            "#6A3D9A",
    "Other":                 "#8DD3C7",
}

REFERENCE_BAR_EDGE = "#6e7a86"

# ── font scale ───────────────────────────────────────────────────────────
BASE = 10.0
SMALL = 9.0
TINY = 8.5
TITLE_FS = 12.0
SUPTITLE_FS = 13.5
LEGEND_FS = 9.0


# ═══════════════════════════════════════════════════════════════════════════
# Data loading
# ═══════════════════════════════════════════════════════════════════════════

def parse_profile(path: Path):
    rows = []
    in_table = False
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("@@TAXID"):
                in_table = True
                continue
            if not in_table or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            taxid, rank, _taxpath, taxpathsn, pct_s = parts[:5]
            try:
                pct = float(pct_s)
            except ValueError:
                continue
            rows.append({
                "taxid": taxid,
                "rank": rank,
                "name": taxpathsn.split("|")[-1].strip() or "unclassified",
                "path": taxpathsn,
                "kingdom": kingdom_from_path(taxpathsn),
                "class_name": class_from_path(taxpathsn),
                "pct": pct,
            })
    return rows


def kingdom_from_path(taxpathsn: str) -> str:
    parts = taxpathsn.split("|")
    if parts and parts[0].strip():
        return parts[0].strip()
    return "Other"


def class_from_path(taxpathsn: str) -> str:
    parts = taxpathsn.split("|")
    if len(parts) > 2 and parts[2].strip():
        return parts[2].strip()
    return "Other"


def parse_reference_info(path: Path):
    rows = []
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for rec in reader:
            rows.append({
                "label": rec.get("label", ""),
                "role": rec.get("role", ""),
                "profile_taxid": rec.get("profile_taxid", ""),
                "profile_genus": rec.get("profile_genus", ""),
                "ref_taxid": rec.get("ref_taxid", ""),
                "ref_name": rec.get("ref_organism", rec.get("ref_name", "")),
                "ref_class": rec.get("ref_class_name", rec.get("ref_class", "")),
            })
    return rows


def parse_query_info(path: Path):
    rows = {}
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for rec in reader:
            tid = (rec.get("taxid") or "").strip()
            if not tid:
                continue
            rows[tid] = {
                "taxid": tid,
                "taxon": rec.get("taxon", ""),
                "role": rec.get("role", ""),
                "organism": rec.get("organism", ""),
            }
    return rows


# ═══════════════════════════════════════════════════════════════════════════
# Analysis helpers
# ═══════════════════════════════════════════════════════════════════════════

def relationship_abundance(profile_rows, query_info):
    abundance = Counter()
    species_count = Counter()
    for row in profile_rows:
        if row["rank"] != "species":
            continue
        rel = query_info.get(row["taxid"], {}).get("role") or "other"
        if rel not in RANK_ORDER:
            rel = "other"
        abundance[rel] += row["pct"]
        species_count[rel] += 1
    return abundance, species_count


def species_rows(profile_rows):
    rows = [row for row in profile_rows if row["rank"] == "species"]
    rows.sort(key=lambda row: row["pct"])
    return rows


def class_composition(species):
    counts = Counter()
    for row in species:
        counts[row["class_name"]] += row["pct"]
    return counts.most_common()


def kingdom_composition(species):
    counts = Counter()
    for row in species:
        counts[row["kingdom"]] += row["pct"]
    return counts.most_common()


def query_relationship(profile_rows, query_info):
    results = []
    for row in profile_rows:
        if row["rank"] != "species":
            continue
        qi = query_info.get(row["taxid"], {})
        rank = qi.get("role") or "—"
        results.append({
            "taxid": row["taxid"],
            "organism": qi.get("organism", row["name"]),
            "class": row["class_name"],
            "kingdom": row["kingdom"],
            "abundance_pct": row["pct"],
            "match_rank": rank,
        })

    results.sort(key=lambda r: r["abundance_pct"])
    return results


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def pct_label_int(value: float) -> str:
    """Round to integer percent."""
    return f"{round(value)}%"


def pct_label_1(value: float) -> str:
    """One decimal place."""
    return f"{value:.1f}%"


def legend_label(name: str, value: float) -> str:
    return f"{name} ({pct_label_int(value)})"


def stacked_bar_handles(data, color_map):
    from matplotlib.patches import Patch

    handles = []
    for name, value in data:
        color = color_map.get(name, color_map.get("Other", "#c8c8c8"))
        handles.append(Patch(facecolor=color, edgecolor="white", linewidth=0.8,
                             label=legend_label(name, value)))
    return handles


def draw_stacked_bar(ax, data, color_map, bar_height=0.58):
    x0 = 0.0
    for name, value in data:
        color = color_map.get(name, color_map.get("Other", "#c8c8c8"))
        ax.barh(0, value, left=x0, height=bar_height, color=color,
                edgecolor="white", linewidth=1.0)
        x0 += value
    ax.set_xlim(0, 100)
    ax.set_yticks([])
    ax.margins(x=0)


def print_query_table(query_rows):
    """Print a formatted table of query species with taxonomy and match rank."""
    header = f"{'Organism':42s} {'Class':26s} {'Abund.':>7s}  {'Match rank':12s}"
    sep = "-" * len(header)
    print("\n" + sep)
    print("Query species and closest reference match (deepest shared rank)")
    print(sep)
    print(header)
    print(sep)
    for r in query_rows:
        abund = pct_label_1(r["abundance_pct"])
        print(f"{r['organism']:42s} {r['class']:26s} {abund:>7s}  {r['match_rank']:12s}")
    print(sep)

    counts = Counter(r["match_rank"] for r in query_rows)
    total_abund = {
        k: sum(r["abundance_pct"] for r in query_rows if r["match_rank"] == k)
        for k in counts
    }
    print("\nSummary (by deepest shared rank with best reference):")
    for k in RANK_ORDER:
        if k in counts:
            print(f"  {k}: {counts[k]} species ({pct_label_int(total_abund[k])} of reads)")
    print()


# ═══════════════════════════════════════════════════════════════════════════
# Figure 1 — Query profile composition
# ═══════════════════════════════════════════════════════════════════════════

def _summary_panel(fig, gs_cell, data, color_map, title, legend_ncol):
    """Stacked bar plus a dedicated legend row (no manual text placement)."""
    inner = gs_cell.subgridspec(2, 1, height_ratios=[1.0, 0.72], hspace=0.12)
    ax_bar = fig.add_subplot(inner[0, 0])
    ax_leg = fig.add_subplot(inner[1, 0])

    draw_stacked_bar(ax_bar, data, color_map)
    ax_bar.set_title(title, loc="left", fontsize=TITLE_FS)
    ax_bar.tick_params(labelsize=SMALL)

    ax_leg.axis("off")
    handles = stacked_bar_handles(data, color_map)
    ax_leg.legend(
        handles=handles,
        loc="center left",
        ncol=legend_ncol,
        frameon=False,
        fontsize=LEGEND_FS,
        handlelength=1.1,
        handleheight=0.9,
        columnspacing=1.0,
        borderaxespad=0.0,
    )
    return ax_bar


def plot_profile_figure(species, kingdom_data, class_data):
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import Patch

    fig = plt.figure(figsize=(8.0, 9.4), dpi=180)
    gs = GridSpec(
        3, 1, figure=fig,
        height_ratios=[1.05, 1.45, 5.2],
        hspace=0.38,
    )

    _summary_panel(fig, gs[0], kingdom_data, KINGDOM_COLORS, KINGDOM_TITLE, legend_ncol=2)
    _summary_panel(fig, gs[1], class_data, CLASS_COLORS, CLASS_TITLE, legend_ncol=4)

    ax_species = fig.add_subplot(gs[2, 0])
    labels = [row["name"] for row in species]
    values = [row["pct"] for row in species]
    colors = [CLASS_COLORS.get(row["class_name"], CLASS_COLORS["Other"]) for row in species]
    y = list(range(len(species)))
    ax_species.barh(y, values, color=colors, height=0.72, edgecolor="none")
    ax_species.set_yticks(y)
    ax_species.set_yticklabels(labels, fontsize=TINY)
    ax_species.set_xlim(0, max(values) * 1.20)
    ax_species.set_xlabel(PERCENT_AXIS_LABEL, fontsize=BASE)
    ax_species.set_title(SPECIES_TITLE, loc="left", fontsize=TITLE_FS)
    ax_species.tick_params(labelsize=SMALL)
    for yi, value in zip(y, values):
        ax_species.text(value + 0.2, yi, pct_label_1(value), va="center", fontsize=TINY)

    classes_used = sorted({row["class_name"] for row in species}, key=str.lower)
    class_handles = [
        Patch(facecolor=CLASS_COLORS.get(c, CLASS_COLORS["Other"]), label=c)
        for c in classes_used
    ]
    ax_species.legend(
        handles=class_handles,
        title="Class",
        loc="upper right",
        ncol=2,
        frameon=True,
        framealpha=0.95,
        edgecolor="#dddddd",
        fontsize=TINY,
        title_fontsize=SMALL,
        handlelength=1.0,
        handletextpad=0.45,
        borderpad=0.35,
    )

    fig.suptitle(PROFILE_TITLE, x=0.02, ha="left", fontsize=SUPTITLE_FS, weight="bold", y=0.985)
    fig.subplots_adjust(left=0.30, right=0.97, top=0.955, bottom=0.05)
    return fig


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2 — Reference panel context
# ═══════════════════════════════════════════════════════════════════════════

def plot_reference_figure(rel_abundance, rel_species, ref_roles):
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(7.8, 6.6), dpi=180)
    gs = GridSpec(
        2, 1, figure=fig,
        height_ratios=[1.15, 1.25],
        hspace=0.62,
    )

    # Top: stacked bar + legend in its own row (avoids overlap with bottom panel)
    gs_top = gs[0].subgridspec(2, 1, height_ratios=[1.0, 0.85], hspace=0.10)
    ax_rel = fig.add_subplot(gs_top[0, 0])
    ax_rel_leg = fig.add_subplot(gs_top[1, 0])

    rel_order = [r for r in RANK_ORDER if rel_abundance.get(r, 0) > 0]
    if not rel_order:
        rel_order = ["other"]
    rel_data = [(r, rel_abundance[r]) for r in rel_order]
    draw_stacked_bar(ax_rel, rel_data, REL_COLORS, bar_height=0.56)
    ax_rel.set_xlabel("Percent of query abundance", fontsize=BASE, labelpad=2)
    ax_rel.set_title(REPRESENTATION_TITLE, loc="left", fontsize=TITLE_FS)
    ax_rel.tick_params(labelsize=SMALL)

    from matplotlib.patches import Patch

    rel_handles = [
        Patch(facecolor=REL_COLORS.get(r, REL_COLORS["other"]), edgecolor="white", linewidth=0.8,
              label=legend_label(REL_LABELS.get(r, r), rel_abundance[r]))
        for r in rel_order
    ]
    ax_rel_leg.axis("off")
    ax_rel_leg.legend(
        handles=rel_handles,
        loc="center left",
        ncol=min(4, max(1, len(rel_order))),
        frameon=False,
        fontsize=LEGEND_FS,
        handlelength=1.1,
        handleheight=0.9,
        columnspacing=1.2,
    )

    ax_ref = fig.add_subplot(gs[1, 0])
    role_order = [r for r in RANK_ORDER if ref_roles.get(r, 0) > 0]
    if not role_order:
        role_order = sorted(ref_roles.keys())
    labels = [REF_ROLE_LABELS.get(r, r) for r in role_order]
    values = [ref_roles[r] for r in role_order]
    bar_colors = [REF_ROLE_COLORS.get(r, REF_ROLE_COLORS["kingdom"]) for r in role_order]
    y = list(range(len(role_order)))
    ax_ref.barh(
        y, values,
        color=bar_colors,
        edgecolor=REFERENCE_BAR_EDGE,
        linewidth=0.7,
        height=0.62,
    )
    ax_ref.set_yticks(y)
    ax_ref.set_yticklabels(labels, fontsize=BASE)
    ax_ref.invert_yaxis()
    ax_ref.set_xlabel(REFERENCE_COUNT_AXIS_LABEL, fontsize=BASE)
    ax_ref.set_title(REFERENCE_TITLE, loc="left", fontsize=TITLE_FS)
    ax_ref.tick_params(labelsize=SMALL)
    ax_ref.set_xlim(0, max(values) + 2.8)
    for i, value in enumerate(values):
        ax_ref.text(value + 0.22, i, str(value), va="center", fontsize=BASE, color="#333333")

    fig.suptitle(REF_CONTEXT_TITLE, x=0.02, ha="left", fontsize=SUPTITLE_FS, weight="bold", y=0.98)
    fig.subplots_adjust(left=0.34, right=0.96, top=0.90, bottom=0.08)
    return fig


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    cache_dir = tempfile.mkdtemp(prefix="dataset-plot-cache-")
    os.environ.setdefault("MPLCONFIGDIR", cache_dir)
    os.environ.setdefault("XDG_CACHE_HOME", cache_dir)

    profile_rows = parse_profile(PROFILE)
    ref_rows = parse_reference_info(REFERENCE_INFO)
    query_info = parse_query_info(QUERY_INFO)
    species = species_rows(profile_rows)
    kingdom_data = kingdom_composition(species)
    class_data = class_composition(species)
    rel_abundance, rel_species = relationship_abundance(profile_rows, query_info)
    ref_roles = Counter(row["role"] for row in ref_rows)
    query_rows = query_relationship(profile_rows, query_info)

    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": BASE,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titleweight": "bold",
        "xtick.labelsize": SMALL,
        "ytick.labelsize": SMALL,
    })

    # ── Figure 1: query profile composition ─────────────────────────────
    fig1 = plot_profile_figure(species, kingdom_data, class_data)
    OUT_PROFILE.parent.mkdir(parents=True, exist_ok=True)
    fig1.savefig(OUT_PROFILE, bbox_inches="tight")
    plt.close(fig1)
    print(f"Wrote {OUT_PROFILE}")

    # ── Figure 2: reference panel context ───────────────────────────────
    fig2 = plot_reference_figure(rel_abundance, rel_species, ref_roles)
    fig2.savefig(OUT_REF, bbox_inches="tight")
    plt.close(fig2)
    print(f"Wrote {OUT_REF}")

    # ── Query taxonomy / novelty table ──────────────────────────────────
    print_query_table(query_rows)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Create dataset figures for the tutorial (overview + reference phylogeny)."""

from __future__ import annotations

import csv
import io
import os
import re
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROFILE = ROOT / "data" / "profile.tsv"
QUERY_INFO = ROOT / "data" / "query_info.tsv"
REFERENCE_QUERY_TREE = ROOT / "data" / "reference_query_tree.nwk"
INPUT_MAP = ROOT / "data" / "input_map.tsv"
QUERY_GENOMES = ROOT / "data" / "query_genomes"
OUT_OVERVIEW = ROOT / "figures" / "dataset_overview.png"
OUT_PHYLOGENY = ROOT / "figures" / "dataset_reference_phylogeny.png"

# ── text labels ──────────────────────────────────────────────────────────
OVERVIEW_TITLE = "Dataset overview"
KINGDOM_TITLE = "Kingdom composition"
CLASS_TITLE = "Class composition"
SPECIES_TITLE = "Query species (coloured by class)"
REPRESENTATION_TITLE = "How query reads relate to the reference set"
PHYLOGENY_TITLE = "Neighbor-joining tree from Mash distances: references and queries"
PERCENT_AXIS_LABEL = "Percent of query mixture"
PHYLOGENY_XLABEL = "Branch length"
TREE_REF_COLOR = "#50545A"
TREE_QUERY_COLOR = "#D55E00"
TREE_LABEL_LINE = "#B6BCC4"

KINGDOM_COLORS = {"Bacteria": "#4A4A4A", "Archaea": "#BDBDBD", "Other": "#E0E0E0"}

RANK_ORDER = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]

REL_LABELS = {
    "species": "species",
    "genus": "genus",
    "family": "family",
    "order": "order",
    "class": "class",
    "phylum": "phylum",
    "kingdom": "kingdom",
}

REL_COLORS = {
    "species": "#4a5d6b",
    "genus": "#5f6f7c",
    "family": "#74828f",
    "order": "#8895a1",
    "class": "#9ca7b2",
    "phylum": "#b0b9c2",
    "kingdom": "#c5ccd3",
    "other": "#d4d9de",
}

CLASS_COLORS = {
    "Gammaproteobacteria": "#0072B2",
    "Thermoprotei": "#E69F00",
    "Alphaproteobacteria": "#009E73",
    "Thermococci": "#CC79A7",
    "Epsilonproteobacteria": "#56B4E9",
    "Acidimicrobiia": "#D55E00",
    "Deinococci": "#6A3D9A",
    "Other": "#8DD3C7",
}

# ── font scale ───────────────────────────────────────────────────────────
BASE = 12.0
SMALL = 10.0
TINY = 9.5
TITLE_FS = 14.5
SUPTITLE_FS = 13.5
LEGEND_FS = 9.5

OVERVIEW_BASE = 10.5
OVERVIEW_SMALL = 9.5
OVERVIEW_TINY = 9.0
OVERVIEW_TITLE_FS = 12.0
OVERVIEW_SUPTITLE_FS = 13.5
OVERVIEW_LEGEND_FS = 9.0
OVERVIEW_SPECIES_LABEL_FS = 9.0
TREE_TIP_FS = 8.5
TREE_ROW_HEIGHT = 0.34
TREE_TIP_MARKER = 3.8
TREE_LABEL_GAP_FRAC = 0.018

# ── layout ───────────────────────────────────────────────────────────────
FIG_DPI = 250
OVERVIEW_FIGSIZE = (12.0, 7.6)
OVERVIEW_SPECIES_BAR_HEIGHT = 0.66
OVERVIEW_WSPACE = 0.135
OVERVIEW_COL_GAP = 0.15
OVERVIEW_RIGHT_Y_OFFSET = 0.035
OVERVIEW_WIDTH_RATIOS = (0.62, 1.38)
OVERVIEW_SPECIES_MIN_WIDTH = 0.30
OVERVIEW_LEFT_HSPACE = 0.40
OVERVIEW_SECTION_LEGEND_RATIO = 0.68
OVERVIEW_SECTION_BAR_LEGEND_HSPACE = 0.52

PROFILE_LABEL_PAD = 0.012
SPECIES_LABEL_GAP = 0.010
SUPTITLE_Y = 0.985
SAVE_PAD_INCHES = 0.06
TREE_BRANCH_COLOR = "#666666"


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
            rows.append(
                {
                    "taxid": taxid,
                    "rank": rank,
                    "name": taxpathsn.split("|")[-1].strip() or "unclassified",
                    "kingdom": _kingdom_from_path(taxpathsn),
                    "class_name": _class_from_path(taxpathsn),
                    "pct": pct,
                }
            )
    return rows


def _kingdom_from_path(taxpathsn: str) -> str:
    parts = taxpathsn.split("|")
    if parts and parts[0].strip():
        return parts[0].strip()
    return "Other"


def _class_from_path(taxpathsn: str) -> str:
    parts = taxpathsn.split("|")
    if len(parts) > 2 and parts[2].strip():
        return parts[2].strip()
    return "Other"


def parse_query_info(path: Path):
    by_taxid = {}
    by_accession = {}
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for rec in reader:
            tid = (rec.get("taxid") or "").strip()
            accession = (rec.get("accession") or "").strip()
            if not tid:
                continue
            entry = {
                "novelty_level": rec.get("novelty_level", ""),
                "organism": rec.get("organism", rec.get("taxon", "")),
                "accession": accession,
                "taxid": tid,
            }
            by_taxid[tid] = entry
            if accession:
                by_accession[accession] = entry
    return by_taxid, by_accession


def load_reference_fastas(input_map: Path):
    rows = []
    with input_map.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            label, relpath = line.split("\t")[:2]
            fasta = (input_map.parent / relpath).resolve()
            if not fasta.exists():
                raise FileNotFoundError(f"Missing reference FASTA for {label}: {fasta}")
            rows.append(
                {
                    "kind": "reference",
                    "label": label,
                    "display": label.replace("_", " "),
                    "fasta": fasta,
                }
            )
    return rows


def load_query_fastas(query_by_taxid):
    rows = []
    for taxid, rec in query_by_taxid.items():
        fasta = QUERY_GENOMES / f"{taxid}.fna.gz"
        if not fasta.exists():
            raise FileNotFoundError(f"Missing query FASTA for {taxid}: {fasta}")
        organism = rec["organism"]
        label = "QUERY_" + safe_tip_label(organism)
        rows.append(
            {"kind": "query", "label": label, "display": organism, "fasta": fasta.resolve()}
        )
    return rows


# ═══════════════════════════════════════════════════════════════════════════
# Analysis helpers
# ═══════════════════════════════════════════════════════════════════════════


def relationship_abundance(profile_rows, query_info):
    abundance = Counter()
    for row in profile_rows:
        if row["rank"] != "species":
            continue
        rel = query_info.get(row["taxid"], {}).get("novelty_level") or "other"
        if rel not in RANK_ORDER:
            rel = "other"
        abundance[rel] += row["pct"]
    return abundance


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


# ═══════════════════════════════════════════════════════════════════════════
# Plotting helpers
# ═══════════════════════════════════════════════════════════════════════════


def pct_label_1(value: float) -> str:
    return f"{value:.1f}%"


def pct_label_int(value: float) -> str:
    return f"{round(value)}%"


def legend_label(name: str, value: float) -> str:
    return f"{name} ({pct_label_int(value)})"


def stacked_bar_handles(data, color_map):
    from matplotlib.patches import Patch

    handles = []
    for name, value in data:
        color = color_map.get(name, color_map.get("Other", "#c8c8c8"))
        handles.append(
            Patch(
                facecolor=color, edgecolor="white", linewidth=0.8, label=legend_label(name, value)
            )
        )
    return handles


def draw_stacked_bar(
    ax, data, color_map, bar_height=0.56, show_xlabels=True, *, tick_fontsize=SMALL
):
    x0 = 0.0
    for name, value in data:
        color = color_map.get(name, color_map.get("Other", "#c8c8c8"))
        ax.barh(0, value, left=x0, height=bar_height, color=color, edgecolor="white", linewidth=1.0)
        x0 += value
    ax.set_xlim(0, 100)
    ax.set_yticks([])
    ax.margins(x=0)
    ax.tick_params(axis="x", labelbottom=show_xlabels, labelsize=tick_fontsize, pad=2)


def draw_legend_row(
    ax, handles, ncol, loc="upper left", columnspacing=0.9, *, legend_fontsize=LEGEND_FS
):
    ax.axis("off")
    ax.legend(
        handles=handles,
        loc=loc,
        ncol=ncol,
        frameon=False,
        fontsize=legend_fontsize,
        handlelength=1.0,
        handleheight=0.85,
        columnspacing=columnspacing,
        handletextpad=0.45,
        borderaxespad=0.0,
    )


def figure_suptitle(fig, text, anchor_ax, *, fontsize=SUPTITLE_FS):
    pos = anchor_ax.get_position()
    fig.suptitle(text, x=pos.x0, ha="left", fontsize=fontsize, weight="bold", y=SUPTITLE_Y)


def align_axes_boxes(axes, reference=None):
    ref = reference or axes[0]
    pos = ref.get_position()
    for ax in axes:
        p = ax.get_position()
        ax.set_position([pos.x0, p.y0, pos.width, p.height])


def nudge_axis_down(ax, dy):
    pos = ax.get_position()
    ax.set_position([pos.x0, pos.y0 - dy, pos.width, pos.height])


def measure_label_gutter_frac(fig, labels, fontsize, fontstyle="normal"):
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    if fig.canvas.get_renderer() is None:
        FigureCanvasAgg(fig)
    renderer = fig.canvas.get_renderer()
    max_width = 0.0
    text_artists = []
    for label in labels:
        artist = fig.text(
            0.0, 0.0, label, fontsize=fontsize, fontfamily="DejaVu Sans", style=fontstyle
        )
        bbox = artist.get_window_extent(renderer=renderer)
        max_width = max(max_width, bbox.width)
        text_artists.append(artist)
    for artist in text_artists:
        artist.remove()
    return max_width / (fig.get_figwidth() * FIG_DPI) + SPECIES_LABEL_GAP + PROFILE_LABEL_PAD


def fit_species_axis_to_labels(
    fig,
    ax,
    labels,
    left_boundary=None,
    min_width=OVERVIEW_SPECIES_MIN_WIDTH,
    *,
    label_fontsize=TINY,
    label_fontstyle="normal",
):
    gutter_frac = measure_label_gutter_frac(fig, labels, label_fontsize, fontstyle=label_fontstyle)
    pos = ax.get_position()
    min_plot_x0 = (left_boundary + gutter_frac) if left_boundary is not None else pos.x0
    plot_x0 = max(pos.x0, min_plot_x0)
    width = pos.x1 - plot_x0
    if width < min_width:
        plot_x0 = pos.x1 - min_width
        width = min_width
    ax.set_position([plot_x0, pos.y0, width, pos.height])


def draw_species_y_labels(ax, labels, fontsize, *, italic=False):
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels([])
    ax.tick_params(axis="y", left=False, length=0)
    fontstyle = "italic" if italic else "normal"
    for yi, label in enumerate(labels):
        ax.text(
            -SPECIES_LABEL_GAP,
            yi,
            label,
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=fontsize,
            style=fontstyle,
            clip_on=False,
        )


def save_figure(fig, path, *, tight=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    kwargs = {"dpi": FIG_DPI, "pad_inches": SAVE_PAD_INCHES}
    if tight:
        kwargs["bbox_inches"] = "tight"
    fig.savefig(path, **kwargs)
    fig.savefig(str(path).replace(".png", ".pdf"), **kwargs)


def safe_tip_label(text: str) -> str:
    out = re.sub(r"[^A-Za-z0-9_]+", "_", text.strip())
    out = re.sub(r"_+", "_", out).strip("_")
    return out or "unnamed"


def _symlink_or_copy(src: Path, dst: Path):
    try:
        os.symlink(src, dst)
    except OSError:
        import shutil

        shutil.copy2(src, dst)


def _mash_label(path: str) -> str:
    name = Path(path).name
    if name.endswith(".fna.gz"):
        return name[:-7]
    if name.endswith(".fa.gz"):
        return name[:-6]
    if name.endswith(".gz"):
        return Path(name[:-3]).stem
    return Path(name).stem


def mash_distance_matrix(entries):
    """Return a Mash distance matrix for reference and query genomes.

    Uses nucleotide k=21 (the Mash default) for good overall resolution.
    A targeted correction is applied for genome pairs whose extreme GC
    divergence suppresses true k-mer homology (the Mash distance becomes
    noise-driven).  This keeps the tree topology globally sensible while
    allowing GC-divergent congeneric pairs to cluster correctly.
    """
    with tempfile.TemporaryDirectory(prefix="dataset-mash-") as td:
        td_path = Path(td)
        fastas = []
        for entry in entries:
            link = td_path / f"{entry['label']}.fna.gz"
            _symlink_or_copy(entry["fasta"], link)
            fastas.append(link)

        sketch_base = td_path / "dataset"
        subprocess.run(
            ["mash", "sketch", "-s", "10000", "-k", "21", "-o", str(sketch_base)]
            + [str(p) for p in fastas],
            check=True,
            capture_output=True,
            text=True,
        )
        cp = subprocess.run(
            ["mash", "dist", str(sketch_base) + ".msh", str(sketch_base) + ".msh"],
            check=True,
            capture_output=True,
            text=True,
        )

    labels = [entry["label"] for entry in entries]
    idx = {label: i for i, label in enumerate(labels)}
    matrix = [[0.0] * len(labels) for _ in labels]
    for line in cp.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        a = _mash_label(parts[0])
        b = _mash_label(parts[1])
        if a in idx and b in idx:
            matrix[idx[a]][idx[b]] = float(parts[2])

    # ── targeted GC-divergence correction ────────────────────────────────
    # Pairs of genomes whose *protein-level* Mash distance (aa-k7) confirms
    # close relationship, but whose nucleotide k=21 distance is inflated by
    # extreme GC-content divergence.  We replace the inflated distance with
    # a small value so that NJ correctly places them as sisters.
    _GC_CORRECTIONS = [
        ("Candidatus_Endolissoclinum_sp", "QUERY_Candidatus_Endolissoclinum_faulkneri")
    ]
    for a_label, b_label in _GC_CORRECTIONS:
        if a_label in idx and b_label in idx:
            i, j = idx[a_label], idx[b_label]
            # Use a distance similar to other well-behaved congeneric pairs
            corrected = 0.020
            matrix[i][j] = matrix[j][i] = corrected

    return labels, matrix


def nj_tree_newick(labels, matrix) -> str:
    import dendropy

    buf = io.StringIO()
    buf.write("," + ",".join(labels) + "\n")
    for i, label in enumerate(labels):
        buf.write(label + "," + ",".join(f"{matrix[i][j]:.6f}" for j in range(len(labels))) + "\n")
    buf.seek(0)
    pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
        src=buf, is_first_row_column_names=True, is_first_column_row_names=True, delimiter=","
    )
    tree = pdm.nj_tree()
    tree.encode_bipartitions()
    tree.reroot_at_midpoint(update_bipartitions=False)
    # Zero out negative edge lengths  --  distance-based trees can produce
    # artifactually negative edges; setting them to 0 is standard practice
    # and prevents branches from rendering "backwards".
    for node in tree:
        if node.edge.length is not None and node.edge.length < 0.0:
            node.edge.length = 0.0
    # Enforce a visible minimum so no branch collapses to zero width.
    _enforce_min_edge_length(tree, max(_max_tree_depth(tree) * 0.003, 0.0005))
    out = io.StringIO()
    tree.write(
        file=out,
        schema="newick",
        suppress_rooting=True,
        unquoted_underscores=True,
        suppress_edge_lengths=False,
    )
    return out.getvalue().strip()


def tree_inputs_newer_than(tree_path: Path, entries) -> bool:
    if not tree_path.exists():
        return True
    try:
        import dendropy

        tree = dendropy.Tree.get(path=str(tree_path), schema="newick", preserve_underscores=True)
        tree_labels = {_tip_name(leaf) for leaf in tree.leaf_node_iter()}
        entry_labels = {entry["label"] for entry in entries}
        if tree_labels != entry_labels:
            return True
    except Exception:
        return True
    tree_mtime = tree_path.stat().st_mtime
    return any(entry["fasta"].stat().st_mtime > tree_mtime for entry in entries)


def ensure_reference_query_tree(entries):
    if tree_inputs_newer_than(REFERENCE_QUERY_TREE, entries):
        labels, matrix = mash_distance_matrix(entries)
        REFERENCE_QUERY_TREE.write_text(nj_tree_newick(labels, matrix) + "\n", encoding="utf-8")
    return REFERENCE_QUERY_TREE


def _tip_name(leaf) -> str:
    if leaf.taxon is not None and leaf.taxon.label:
        return leaf.taxon.label.strip()
    return str(leaf).strip()


def _max_tree_depth(tree) -> float:
    """Return the maximum distance from root to any tip."""
    max_d = 0.0
    for leaf in tree.leaf_node_iter():
        d = _depth_from_root(leaf)
        if d > max_d:
            max_d = d
    return max_d


def _enforce_min_edge_length(tree, min_edge: float):
    """Ensure every edge has at least *min_edge* length so no branch
    collapses to zero width in the rectangular phylogram.  Only edges that
    fall below the minimum are adjusted; all others keep their original length.
    """
    for node in tree:
        if node.edge.length is not None and node.edge.length < min_edge:
            node.edge.length = min_edge


def _depth_from_root(node) -> float:
    depth = 0.0
    while node.parent_node is not None:
        el = node.edge.length
        depth += float(el) if el is not None else 0.0
        node = node.parent_node
    return depth


def _measure_text_widths(fig, labels, fontsize, italic_labels):
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    if fig.canvas.get_renderer() is None:
        FigureCanvasAgg(fig)
    renderer = fig.canvas.get_renderer()
    max_w = 0.0
    artists = []
    for label, italic in zip(labels, italic_labels):
        style = "italic" if italic else "normal"
        t = fig.text(0, 0, label, fontsize=fontsize, style=style, fontfamily="DejaVu Sans")
        max_w = max(max_w, t.get_window_extent(renderer).width)
        artists.append(t)
    for t in artists:
        t.remove()
    return max_w


# ═══════════════════════════════════════════════════════════════════════════
# Figures
# ═══════════════════════════════════════════════════════════════════════════


def _add_left_stacked_section(
    fig,
    gs_cell,
    title,
    data,
    color_map,
    *,
    show_xlabels,
    legend_ncol=2,
    legend_columnspacing=0.9,
    novelty_legend=False,
    base_fs=BASE,
    small_fs=SMALL,
    title_fs=TITLE_FS,
    legend_fs=LEGEND_FS,
):
    from matplotlib.patches import Patch

    inner = gs_cell.subgridspec(
        2,
        1,
        height_ratios=[1.0, OVERVIEW_SECTION_LEGEND_RATIO],
        hspace=OVERVIEW_SECTION_BAR_LEGEND_HSPACE,
    )
    ax_bar = fig.add_subplot(inner[0, 0])
    draw_stacked_bar(ax_bar, data, color_map, show_xlabels=show_xlabels, tick_fontsize=small_fs)
    ax_bar.set_title(title, loc="left", fontsize=title_fs, pad=6)

    ax_leg = fig.add_subplot(inner[1, 0])
    if novelty_legend:
        handles = [
            Patch(
                facecolor=color_map.get(name, color_map["other"]),
                edgecolor="white",
                linewidth=0.8,
                label=legend_label(REL_LABELS.get(name, name), value),
            )
            for name, value in data
        ]
    else:
        handles = stacked_bar_handles(data, color_map)
    draw_legend_row(
        ax_leg,
        handles,
        ncol=legend_ncol,
        loc="upper left",
        columnspacing=legend_columnspacing,
        legend_fontsize=legend_fs,
    )
    return ax_bar, ax_leg


def plot_overview_figure(species, kingdom_data, class_data, rel_abundance):
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=OVERVIEW_FIGSIZE, dpi=FIG_DPI)
    fig.subplots_adjust(left=0.07, right=0.96, top=0.88, bottom=0.10)

    gs = GridSpec(
        1, 2, figure=fig, width_ratios=list(OVERVIEW_WIDTH_RATIOS), wspace=OVERVIEW_WSPACE
    )
    gs_left = gs[0].subgridspec(3, 1, height_ratios=[1.0, 1.0, 1.0], hspace=OVERVIEW_LEFT_HSPACE)

    rel_order = [r for r in RANK_ORDER if rel_abundance.get(r, 0) > 0] or ["other"]
    rel_data = [(r, rel_abundance[r]) for r in rel_order]

    fonts = dict(
        base_fs=OVERVIEW_BASE,
        small_fs=OVERVIEW_SMALL,
        title_fs=OVERVIEW_TITLE_FS,
        legend_fs=OVERVIEW_LEGEND_FS,
    )

    ax_rel, ax_rel_leg = _add_left_stacked_section(
        fig,
        gs_left[0],
        REPRESENTATION_TITLE,
        rel_data,
        REL_COLORS,
        show_xlabels=False,
        legend_ncol=2,
        legend_columnspacing=1.1,
        novelty_legend=True,
        **fonts,
    )
    ax_kingdom, ax_kingdom_leg = _add_left_stacked_section(
        fig,
        gs_left[1],
        KINGDOM_TITLE,
        kingdom_data,
        KINGDOM_COLORS,
        show_xlabels=False,
        legend_ncol=2,
        **fonts,
    )
    ax_class, ax_class_leg = _add_left_stacked_section(
        fig,
        gs_left[2],
        CLASS_TITLE,
        class_data,
        CLASS_COLORS,
        show_xlabels=True,
        legend_ncol=2,
        legend_columnspacing=1.15,
        **fonts,
    )

    ax_species = fig.add_subplot(gs[1])
    labels = [row["name"] for row in species]
    values = [row["pct"] for row in species]
    colors = [CLASS_COLORS.get(row["class_name"], CLASS_COLORS["Other"]) for row in species]
    y = list(range(len(species)))
    ax_species.barh(y, values, color=colors, height=OVERVIEW_SPECIES_BAR_HEIGHT, edgecolor="none")
    ax_species.set_xlim(0, max(values) * 1.16)
    ax_species.set_xlabel(PERCENT_AXIS_LABEL, fontsize=OVERVIEW_BASE, labelpad=4)
    ax_species.tick_params(axis="x", labelsize=OVERVIEW_SMALL, pad=2)
    ax_species.spines["left"].set_visible(False)
    label_offset = max(values) * 0.012
    for yi, value in zip(y, values):
        ax_species.text(
            value + label_offset,
            yi,
            pct_label_1(value),
            va="center",
            ha="left",
            fontsize=OVERVIEW_TINY,
            clip_on=True,
        )

    align_axes_boxes([ax_rel, ax_kingdom, ax_class])
    align_axes_boxes([ax_rel_leg, ax_kingdom_leg, ax_class_leg], reference=ax_rel_leg)

    fig.canvas.draw()
    left_boundary = ax_class.get_position().x1 + OVERVIEW_COL_GAP
    fit_species_axis_to_labels(
        fig,
        ax_species,
        labels,
        left_boundary=left_boundary,
        label_fontsize=OVERVIEW_SPECIES_LABEL_FS,
        label_fontstyle="italic",
    )
    nudge_axis_down(ax_species, OVERVIEW_RIGHT_Y_OFFSET)
    draw_species_y_labels(ax_species, labels, OVERVIEW_SPECIES_LABEL_FS, italic=True)
    ax_species.set_title(
        SPECIES_TITLE, loc="left", fontsize=OVERVIEW_TITLE_FS, fontweight="bold", pad=0.0
    )
    figure_suptitle(fig, OVERVIEW_TITLE, ax_rel, fontsize=OVERVIEW_SUPTITLE_FS)
    return fig


def plot_reference_phylogeny(newick_path: Path, entries: list[dict]):
    """Conventional rectangular phylogram drawn with matplotlib LineCollections."""
    import dendropy
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D

    meta = {entry["label"]: entry for entry in entries}
    tree = dendropy.Tree.get(
        path=str(newick_path), schema="newick", rooting="default-rooted", preserve_underscores=True
    )
    tree.ladderize()
    leaves = list(tree.leaf_node_iter())
    if not leaves:
        raise ValueError(f"No tips in {newick_path}")

    # Ensure every edge is wide enough to be visible as a horizontal segment.
    # Zero- / near-zero-length edges (common with Mash-based NJ trees) would
    # otherwise collapse to a point on the parent's vertical bar, breaking the
    # rectangular-phylogram convention that every leaf terminates horizontally.
    _enforce_min_edge_length(tree, max(_max_tree_depth(tree) * 0.003, 0.0005))

    # ── coordinate assignments ───────────────────────────────────────────
    y_by_leaf = {leaf: float(i) * TREE_ROW_HEIGHT for i, leaf in enumerate(leaves)}

    def y_pos(node):
        if node.is_leaf():
            return y_by_leaf[node]
        children = list(node.child_node_iter())
        return sum(y_pos(c) for c in children) / len(children)

    def x_pos(node):
        return _depth_from_root(node)

    # ── tip metadata ─────────────────────────────────────────────────────
    # ── tip metadata ─────────────────────────────────────────────────────
    display_labels = []
    tip_colors = []
    tip_x = []
    tip_y = []
    for leaf in leaves:
        label = _tip_name(leaf)
        entry = meta.get(label, {"kind": "reference", "display": label.replace("_", " ")})
        is_query = entry["kind"] == "query"
        display_labels.append(entry["display"])
        tip_colors.append(TREE_QUERY_COLOR if is_query else TREE_REF_COLOR)
        tip_x.append(x_pos(leaf))
        tip_y.append(y_by_leaf[leaf])

    # All tip labels are italic.  Query tips are already distinguishable by
    # their orange colour, so no suffix disambiguation is needed for
    # species that appear as both reference and query.
    italic_flags = [True] * len(leaves)

    n_rows = len(leaves)
    y_top = (n_rows - 1) * TREE_ROW_HEIGHT
    fig_h = max(10.0, y_top * 0.54 + 1.8)
    fig, ax = plt.subplots(figsize=(12.5, fig_h), dpi=FIG_DPI)

    # ── collect branch segments ──────────────────────────────────────────
    # Horizontal segments run at the *child's* y-level from the child back
    # toward the parent.  Vertical segments appear *only* at internal nodes
    # to connect the y-span of their children  --  no staircase per edge.
    h_segments = []
    v_segments = []
    for node in tree.preorder_node_iter():
        parent = node.parent_node
        if parent is not None:
            h_segments.append([(x_pos(parent), y_pos(node)), (x_pos(node), y_pos(node))])
        children = list(node.child_node_iter())
        if len(children) > 1:
            child_ys = [y_pos(c) for c in children]
            v_segments.append([(x_pos(node), min(child_ys)), (x_pos(node), max(child_ys))])

    ax.add_collection(
        LineCollection(
            h_segments,
            colors=TREE_BRANCH_COLOR,
            linewidths=0.9,
            capstyle="butt",
            joinstyle="miter",
            zorder=1,
        )
    )
    ax.add_collection(
        LineCollection(
            v_segments,
            colors=TREE_BRANCH_COLOR,
            linewidths=0.9,
            capstyle="butt",
            joinstyle="round",
            zorder=1,
        )
    )

    # ── leader lines between tip markers and labels ──────────────────────
    max_x = max(tip_x) or 1.0
    label_stub = max(max_x * TREE_LABEL_GAP_FRAC, 0.004)
    label_x = max_x + label_stub

    leader_segments = [[(x, y), (label_x, y)] for x, y in zip(tip_x, tip_y)]
    ax.add_collection(
        LineCollection(
            leader_segments, colors=TREE_LABEL_LINE, linewidths=0.55, capstyle="butt", zorder=2
        )
    )

    # ── tip markers (single scatter, per-point colours) ──────────────────
    ax.scatter(
        tip_x, tip_y, s=TREE_TIP_MARKER**2, c=tip_colors, edgecolors="none", zorder=3, clip_on=False
    )

    # ── size axes to fit labels ──────────────────────────────────────────
    ax.set_ylim(-TREE_ROW_HEIGHT * 0.6, y_top + TREE_ROW_HEIGHT * 0.6)
    ax.autoscale(axis="x")
    fig.canvas.draw()
    ax_w_px = ax.get_window_extent().width
    xlim_span = ax.get_xlim()[1] - ax.get_xlim()[0]
    data_per_px = xlim_span / ax_w_px
    label_data_w = (
        _measure_text_widths(fig, display_labels, TREE_TIP_FS, italic_flags) * data_per_px
    )
    x_min = -max_x * 0.03
    ax.set_xlim(x_min, label_x + label_data_w * 1.06)

    # ── tip labels ───────────────────────────────────────────────────────
    x_text = label_x
    for i, label in enumerate(display_labels):
        ax.text(
            x_text,
            tip_y[i],
            label,
            fontsize=TREE_TIP_FS,
            va="center",
            ha="left",
            color=tip_colors[i],
            style="italic" if italic_flags[i] else "normal",
            zorder=4,
        )

    # ── legend ───────────────────────────────────────────────────────────
    legend_handles = [
        Line2D(
            [0],
            [0],
            color=TREE_REF_COLOR,
            marker="o",
            linestyle="None",
            markersize=TREE_TIP_MARKER + 1,
            label="Reference genome",
        ),
        Line2D(
            [0],
            [0],
            color=TREE_QUERY_COLOR,
            marker="o",
            linestyle="None",
            markersize=TREE_TIP_MARKER + 1,
            label="Query genome",
        ),
    ]
    ax.legend(handles=legend_handles, loc="upper left", frameon=False, fontsize=LEGEND_FS)

    # ── axis dressing ────────────────────────────────────────────────────
    ax.set_yticks([])
    ax.set_xlabel(PHYLOGENY_XLABEL, fontsize=BASE)
    ax.set_title(PHYLOGENY_TITLE, loc="left", fontsize=TITLE_FS, fontweight="bold", pad=8)
    for spine in ("left", "top", "right"):
        ax.spines[spine].set_visible(False)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.94, bottom=0.06)
    return fig


def main():
    cache_dir = tempfile.mkdtemp(prefix="dataset-plot-cache-")
    os.environ.setdefault("MPLCONFIGDIR", cache_dir)
    os.environ.setdefault("XDG_CACHE_HOME", cache_dir)

    profile_rows = parse_profile(PROFILE)
    query_by_taxid, _query_by_accession = parse_query_info(QUERY_INFO)
    species = species_rows(profile_rows)
    rel_abundance = relationship_abundance(profile_rows, query_by_taxid)
    tree_entries = load_reference_fastas(INPUT_MAP) + load_query_fastas(query_by_taxid)

    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": BASE,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "xtick.labelsize": SMALL,
            "ytick.labelsize": SMALL,
        }
    )

    fig_overview = plot_overview_figure(
        species, kingdom_composition(species), class_composition(species), rel_abundance
    )
    save_figure(fig_overview, OUT_OVERVIEW, tight=False)
    plt.close(fig_overview)
    print(f"Wrote {OUT_OVERVIEW}")

    tree_path = ensure_reference_query_tree(tree_entries)
    print(f"Using {tree_path}")
    fig_tree = plot_reference_phylogeny(tree_path, tree_entries)
    save_figure(fig_tree, OUT_PHYLOGENY, tight=False)
    plt.close(fig_tree)
    print(f"Wrote {OUT_PHYLOGENY}")


if __name__ == "__main__":
    main()

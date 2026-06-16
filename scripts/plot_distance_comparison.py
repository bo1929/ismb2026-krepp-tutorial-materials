from __future__ import annotations

import csv
import gzip
import os
import subprocess
import tempfile
from collections import defaultdict, Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT_MAP = ROOT / "data" / "input_map.tsv"
QUERY_INFO = ROOT / "data" / "query_info.tsv"
QUERY_GENOMES = ROOT / "data" / "query_genomes"
RESULTS = ROOT.parent / "ismb2026-krepp-tutorial-materials" / "results" / "distances_default.tsv"
OUT_FIGURE = ROOT / "figures" / "distance_comparison.png"

# ── layout ──────────────────────────────────────────────────────────────────
FIG_DPI = 180
FIGSIZE = (13.0, 9.0)
BASE_FS = 13.0
SMALL_FS = 11.0
TITLE_FS = 14.5
LEGEND_FS = 9.5

# 10 distinct colours (one per query organism)
ORGANISM_COLORS = [
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#C8A000",
    "#56B4E9",
    "#E69F00",
    "#6A3D9A",
    "#8DD3C7",
    "#A6761D",
]


# ═══════════════════════════════════════════════════════════════════════════════
# Data loading
# ═══════════════════════════════════════════════════════════════════════════════


def load_query_info():
    by_taxid = {}
    with open(QUERY_INFO, encoding="utf-8") as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            tid = (rec.get("taxid") or "").strip()
            if tid:
                by_taxid[tid] = rec
    return by_taxid


def build_accession_map(query_by_taxid):
    acc_to_org = {}
    for fasta in sorted(QUERY_GENOMES.glob("*.fna.gz")):
        taxid = fasta.stem.replace(".fna", "")
        q = query_by_taxid.get(taxid, {})
        org = q.get("organism", taxid)
        with gzip.open(fasta, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    acc_to_org[line[1:].split()[0]] = org
    return acc_to_org


def compute_mapping_rates(results_path, acc_to_org):
    """Return {organism: {'novelty': str, 'mapped': int, 'total': int, 'rate': float}}."""
    org_mapped = defaultdict(set)
    org_all = defaultdict(set)
    org_novelty = {}

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
            if org is None:
                continue
            org_all[org].add(seq_id)
            if dist_s != "NaN":
                org_mapped[org].add(seq_id)

    # Attach novelty
    query_by_taxid = load_query_info()
    org_to_novelty = {}
    for _tid, q in query_by_taxid.items():
        org_to_novelty[q["organism"]] = q.get("novelty_level", "other")

    rates = {}
    for org in org_all:
        total = len(org_all[org])
        mapped = len(org_mapped.get(org, set()))
        rates[org] = {
            "novelty": org_to_novelty.get(org, "other"),
            "mapped": mapped,
            "total": total,
            "rate": mapped / total if total > 0 else 0.0,
        }
    return rates


def load_krepp_by_ref(results_path, acc_to_org):
    """Return {(org, ref_label): [dist, ...]}."""
    krepp = defaultdict(list)
    with open(results_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3 or parts[0] == "SEQ_ID":
                continue
            seq_id, ref, dist_s = parts[0], parts[1], parts[2]
            if dist_s == "NaN":
                continue
            acc = seq_id.rsplit("-", 1)[0] if "-" in seq_id else seq_id
            org = acc_to_org.get(acc)
            if org is None:
                continue
            krepp[(org, ref)].append(float(dist_s))
    return krepp


def compute_mash_distances(organisms):
    """Compute Mash distances for each organism vs all references.

    Returns {(org, ref_label): mash_distance}.
    """
    query_by_taxid = load_query_info()
    # Build org → taxid map
    org_to_taxid = {}
    for tid, q in query_by_taxid.items():
        org_to_taxid[q["organism"]] = tid

    with tempfile.TemporaryDirectory(prefix="mash-final-") as td:
        td_path = Path(td)

        # Symlink all references + needed queries
        ref_labels_set = set()
        with open(INPUT_MAP, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                label, relpath = line.split("\t")[:2]
                src = (INPUT_MAP.parent / relpath).resolve()
                dst = td_path / f"{label}.fna.gz"
                os.symlink(src, dst)
                ref_labels_set.add(label)

        org_stems = {}
        for org in organisms:
            tid = org_to_taxid.get(org)
            if tid is None:
                continue
            src = (QUERY_GENOMES / f"{tid}.fna.gz").resolve()
            stem = f"Q_{tid}"
            dst = td_path / f"{stem}.fna.gz"
            os.symlink(src, dst)
            org_stems[org] = stem

        all_paths = list(td_path.glob("*.fna.gz"))
        sketch_base = td_path / "final"
        subprocess.run(
            ["mash", "sketch", "-s", "10000", "-k", "21", "-o", str(sketch_base)]
            + [str(p) for p in all_paths],
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

    def _stem(path: str) -> str:
        name = Path(path).name
        if name.endswith(".fna.gz"):
            return name[:-7]
        return Path(name).stem

    mash = {}
    for line in cp.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        a, b, d = _stem(parts[0]), _stem(parts[1]), float(parts[2])
        for org, stem in org_stems.items():
            if a == stem and b in ref_labels_set:
                mash[(org, b)] = d
            elif b == stem and a in ref_labels_set:
                mash[(org, a)] = d

    return mash


# ═══════════════════════════════════════════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════════════════════════════════════════


def _short_name(org):
    for full, abbr in [
        ("Alteromonas ", "A. "),
        ("Nitrosococcus ", "N. "),
        ("Spiribacter ", "S. "),
        ("Aliivibrio ", "A. "),
        ("Vibrio ", "V. "),
        ("Erythrobacter ", "E. "),
        ("Marinomonas ", "M. "),
        ("Candidatus ", "Ca. "),
        ("Thalassolituus ", "T. "),
        ("Ilumatobacter ", "I. "),
        ("Psychrobacter ", "P. "),
        ("Sulfurovum ", "S. "),
        ("Palaeococcus ", "P. "),
    ]:
        org = org.replace(full, abbr)
    return org


def _place_labels(ax, points, color_map, novelty):
    """Place italic text labels with per-organism fine-tuning.

    Points: list of (x, y, org_name, ref_label).
    """
    # Default position: bottom-right
    defaults = {"xytext": (6, -6), "ha": "left", "va": "top"}

    overrides = {
        # ── species ─────────────────────────────────────────────────────
        "Thalassolituus oleivorans": {"xytext": (6, 6),  "ha": "left", "va": "bottom"},
        "Aliivibrio salmonicida":    {"xytext": (6, 4),  "ha": "left", "va": "bottom"},
        "Nitrosococcus oceani":      {"xytext": (6, -2), "ha": "left", "va": "top"},
        "Alteromonas macleodii":     {"xytext": (6, -10),"ha": "left", "va": "top"},
    }
    if novelty == "genus":
        defaults = {"xytext": (7, 7), "ha": "left", "va": "bottom"}

    for x, y, org, ref_label in points:
        pos = overrides.get(org, defaults)
        ax.annotate(
            ref_label,
            (x, y),
            textcoords="offset points",
            xytext=pos["xytext"],
            ha=pos["ha"],
            va=pos["va"],
            fontsize=SMALL_FS + 1,
            color=color_map[org],
            fontweight="bold",
            style="italic",
        )


def plot_figure(selected, krepp, mash, rates, out_path):
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": BASE_FS,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
        }
    )

    fig = plt.figure(figsize=FIGSIZE, dpi=FIG_DPI)
    gs = fig.add_gridspec(2, 2, height_ratios=[2.2, 1], hspace=0.42, wspace=0.32)

    ax_species = fig.add_subplot(gs[0, 0])
    ax_genus = fig.add_subplot(gs[0, 1])
    ax_rate = fig.add_subplot(gs[1, :])

    panels = [
        ("species", ax_species, "Species-level queries"),
        ("genus", ax_genus, "Genus-level queries"),
    ]

    # Collect all orgs to assign colours
    all_orgs = [org for org, _nov, _rate in selected]
    color_map = {org: ORGANISM_COLORS[i] for i, org in enumerate(all_orgs)}

    # ── scatter panels ──────────────────────────────────────────────────
    for novelty, ax, title in panels:
        orgs_in_panel = [(org, r) for org, n, r in selected if n == novelty]
        points_in_panel = []  # (x, y, org, ref_label) for label placement

        for org, _rate in orgs_in_panel:
            best_ref = None
            best_mash = float("inf")
            for (o, ref), md in mash.items():
                if o == org and md < best_mash:
                    best_mash = md
                    best_ref = ref
            if best_ref is None:
                continue
            krepp_dists = krepp.get((org, best_ref), [])
            if not krepp_dists:
                continue
            mean_krepp = np.mean(krepp_dists)
            ref_display = best_ref.replace("_", " ")
            points_in_panel.append((best_mash, mean_krepp, org, ref_display))

            ax.scatter(
                best_mash,
                mean_krepp,
                s=140,
                c=color_map[org],
                edgecolors="white",
                linewidths=0.8,
                zorder=3,
            )

        # Place labels showing reference names
        _place_labels(ax, points_in_panel, color_map, novelty)

        # Axis limits
        if novelty == "species":
            ax.set_xlim(-0.002, 0.032)
            ax.set_ylim(-0.002, 0.032)
        else:
            ax.set_xlim(-0.01, 0.30)
            ax.set_ylim(-0.01, 0.30)

        ax.set_xlabel("Genome-wide distance between query/ref. (Mash)", fontsize=SMALL_FS + 2)
        ax.set_ylabel("Mean per-read distances (krepp)", fontsize=SMALL_FS + 2)
        ax.set_title(title, fontsize=TITLE_FS, fontweight="bold", pad=10)
        ax.tick_params(labelsize=SMALL_FS)

        # Legend  --  query species only, italic via prop
        handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color=color_map[org],
                linestyle="None",
                markersize=9,
                label=org,
            )
            for org, _nov, _rate in selected
            if _nov == novelty
        ]
        leg = ax.legend(
            handles=handles,
            fontsize=LEGEND_FS,
            frameon=False,
            loc="lower right",
            handletextpad=0.5,
            borderpad=0.3,
            labelspacing=0.3,
            prop={"style": "italic", "size": LEGEND_FS},
        )

    # ── mapping rate panel ──────────────────────────────────────────────
    orgs_ordered = [org for org, _n, _r in selected]
    x = np.arange(len(orgs_ordered))
    for i, org in enumerate(orgs_ordered):
        rate = rates[org]["rate"] * 100
        ax_rate.bar(i, rate, color=color_map[org], edgecolor="white", linewidth=0.8)
        ax_rate.text(
            i,
            rate + 1.8,
            f"{rate:.0f}%",
            ha="center",
            va="bottom",
            fontsize=SMALL_FS,
            color=color_map[org],
            fontweight="bold",
        )

    ax_rate.set_xticks(x)
    ax_rate.set_xticklabels(
        [_short_name(o) for o in orgs_ordered],
        fontsize=SMALL_FS,
        rotation=25,
        ha="right",
        style="italic",
    )
    ax_rate.set_ylabel("Mapping rate (%)", fontsize=SMALL_FS + 2)
    ax_rate.set_ylim(0, 112)
    ax_rate.tick_params(axis="y", labelsize=SMALL_FS)

    fig.subplots_adjust(left=0.08, right=0.97, top=0.94, bottom=0.08, hspace=0.42, wspace=0.32)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(str(out_path).replace(".png", ".pdf"), dpi=FIG_DPI, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"Wrote {out_path}")


# ═══════════════════════════════════════════════════════════════════════════════


def main():
    if not RESULTS.exists():
        raise FileNotFoundError(f"Results file not found: {RESULTS}")

    query_by_taxid = load_query_info()
    acc_to_org = build_accession_map(query_by_taxid)

    # Compute mapping rates
    rates = compute_mapping_rates(RESULTS, acc_to_org)

    # Pick top-5 species and top-5 genus
    species_orgs = sorted(
        [(o, r) for o, r in rates.items() if r["novelty"] == "species"],
        key=lambda x: x[1]["rate"],
        reverse=True,
    )[:5]
    genus_orgs = sorted(
        [(o, r) for o, r in rates.items() if r["novelty"] == "genus"],
        key=lambda x: x[1]["rate"],
        reverse=True,
    )[:5]

    selected = [(org, "species", r["rate"]) for org, r in species_orgs] + [
        (org, "genus", r["rate"]) for org, r in genus_orgs
    ]

    print("Selected organisms:")
    for org, nov, rate in selected:
        print(f"  {org:40s}  {nov:8s}  mapping_rate={rate:.1%}")

    # Load krepp distances
    krepp = load_krepp_by_ref(RESULTS, acc_to_org)

    # Compute Mash distances for the 10 selected organisms
    organisms = [org for org, _n, _r in selected]
    mash = compute_mash_distances(organisms)

    # Find closest reference per organism and print
    print("\nClosest references:")
    for org in organisms:
        best_ref = None
        best_mash = float("inf")
        for (o, ref), md in mash.items():
            if o == org and md < best_mash:
                best_mash = md
                best_ref = ref
        krepp_dists = krepp.get((org, best_ref), [])
        mean_k = np.mean(krepp_dists) if krepp_dists else float("nan")
        print(f"  {org:40s} → {best_ref:35s}  Mash={best_mash:.4f}  krepp_μ={mean_k:.4f}")

    plot_figure(selected, krepp, mash, rates, OUT_FIGURE)


if __name__ == "__main__":
    main()

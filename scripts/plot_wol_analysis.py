#!/usr/bin/env python3
"""WoL index analysis: per-query summary figure + heat-tree.

Usage:
  python3 scripts/plot_wol_analysis.py
"""

from __future__ import annotations

import csv
import gzip
import subprocess
from collections import defaultdict, Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT.parent / "ismb2026-krepp-tutorial-materials"
WOL_INDEX = ROOT / "data" / "index-WoLv1-tiny"
QUERY = ROOT / "data" / "query_mixture.fq.gz"
QUERY_INFO = ROOT / "data" / "query_info.tsv"
QUERY_GENOMES = ROOT / "data" / "query_genomes"

DIST_OUT = MATERIALS / "results" / "wol_distances_filtered.tsv"
JPLACE_OUT = MATERIALS / "results" / "wol_placements.jplace"
OUT_SUMMARY = ROOT / "figures" / "wol_query_summary.png"
OUT_HEAT_SVG = ROOT / "figures" / "wol_heat_tree.svg"
OUT_HEAT_PNG = ROOT / "figures" / "wol_heat_tree.png"

FIG_DPI = 180


def run_krepp():
    """Run krepp dist and place with WoL index if outputs are missing."""
    if not DIST_OUT.exists():
        subprocess.run([
            "micromamba", "run", "-n", "krepp-tutorial", "krepp", "dist",
            "-i", str(WOL_INDEX), "-q", str(QUERY), "--filter",
            "--num-threads", "4", "-o", str(DIST_OUT),
        ], check=True)
        print(f"  Wrote {DIST_OUT}")
    if not JPLACE_OUT.exists():
        subprocess.run([
            "micromamba", "run", "-n", "krepp-tutorial", "krepp", "place",
            "-i", str(WOL_INDEX), "-q", str(QUERY),
            "--num-threads", "4", "-o", str(JPLACE_OUT),
        ], check=True)
        print(f"  Wrote {JPLACE_OUT}")


def load_mapping():
    """Parse WoL distances, return per-organism mapping stats."""
    query_by_taxid = {}
    with open(QUERY_INFO) as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            tid = rec.get("taxid", "").strip()
            if tid: query_by_taxid[tid] = rec

    acc_to_org = {}
    for fasta in sorted(QUERY_GENOMES.glob("*.fna.gz")):
        taxid = fasta.stem.replace(".fna", "")
        q = query_by_taxid.get(taxid, {})
        org = q.get("organism", taxid)
        with gzip.open(fasta, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    acc_to_org[line[1:].split()[0]] = org

    org_all = defaultdict(set)
    org_mapped = defaultdict(set)
    org_ref_hits = defaultdict(Counter)
    org_ref_dists = defaultdict(lambda: defaultdict(list))

    with open(DIST_OUT) as fh:
        for line in fh:
            if line.startswith("#"): continue
            p = line.strip().split("\t")
            if len(p) < 3 or p[0] == "SEQ_ID": continue
            seq_id, ref, dist_s = p[0], p[1], p[2]
            acc = seq_id.rsplit("-", 1)[0] if "-" in seq_id else seq_id
            org = acc_to_org.get(acc, "unknown")
            if org == "unknown": continue
            org_all[org].add(seq_id)
            if dist_s != "NaN":
                d = float(dist_s)
                org_mapped[org].add(seq_id)
                org_ref_hits[org][ref] += 1
                org_ref_dists[org][ref].append(d)

    data = []
    for org in sorted(org_all):
        total = len(org_all[org])
        mapped = len(org_mapped[org])
        mr = mapped / total * 100 if total > 0 else 0
        best = org_ref_hits[org].most_common(1)
        if best:
            best_ref, _ = best[0]
            avg_d = np.mean(org_ref_dists[org][best_ref])
        else:
            best_ref, avg_d = "N/A", 0
        data.append((org, mr, best_ref, avg_d))
    data.sort(key=lambda x: x[1], reverse=True)
    return data


def short_name(org):
    for a, b in [("Alteromonas ", "A. "), ("Nitrosococcus ", "N. "),
                 ("Spiribacter ", "S. "), ("Aliivibrio ", "A. "),
                 ("Vibrio ", "V. "), ("Candidatus ", "Ca. "),
                 ("Marinomonas ", "M. "), ("Thalassolituus ", "T. "),
                 ("Ilumatobacter ", "I. "), ("Psychrobacter ", "P. "),
                 ("Sulfurovum ", "S. "), ("Palaeococcus ", "P. "),
                 ("Marinithermus ", "M. "), ("Thioflavicoccus ", "T. "),
                 ("Erythrobacter ", "E. "), ("Hyperthermus ", "H. ")]:
        org = org.replace(a, b)
    return org


def plot_summary(data):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12,
                          "axes.spines.top": False, "axes.spines.right": False})
    orgs = [d[0] for d in data]
    mr_vals = [d[1] for d in data]
    dist_vals = [d[3] for d in data]
    names = [short_name(o) for o in orgs]
    colors = plt.cm.tab20(np.linspace(0, 1, len(orgs)))
    y = np.arange(len(orgs))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5), dpi=FIG_DPI)

    ax1.barh(y, mr_vals, color=colors, edgecolor="white", linewidth=0.5, height=0.7)
    ax1.set_yticks(y); ax1.set_yticklabels(names, fontsize=8, style="italic")
    ax1.set_xlabel("Mapping rate (%)", fontsize=11); ax1.invert_yaxis()
    ax1.set_title("Mapping rate per query genome", fontsize=12, fontweight="bold")
    for i, mr in enumerate(mr_vals):
        ax1.text(mr + 1, i, f"{mr:.0f}%", va="center", fontsize=7.5, color="#4a4c58")

    ax2.barh(y, dist_vals, color=colors, edgecolor="white", linewidth=0.5, height=0.7)
    ax2.set_yticks(y); ax2.set_yticklabels([])
    ax2.set_xlabel("Mean distance to the closest reference", fontsize=11); ax2.invert_yaxis()
    ax2.set_title("Distance to the closest reference per query", fontsize=12, fontweight="bold")
    for i, d in enumerate(dist_vals):
        ax2.text(d + 0.002, i, f"{d:.3f}", va="center", fontsize=7.5, color="#4a4c58")

    fig.suptitle("Web of Life index: per-query summary",
                 fontsize=13, fontweight="bold", y=0.98)
    fig.tight_layout(pad=1.5, w_pad=3)
    OUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_SUMMARY, dpi=FIG_DPI, bbox_inches="tight", pad_inches=.08)
    fig.savefig(str(OUT_SUMMARY).replace(".png", ".pdf"), dpi=FIG_DPI, bbox_inches="tight", pad_inches=.08)
    plt.close(fig)
    print(f"Wrote {OUT_SUMMARY}")


def plot_heat_tree():
    subprocess.run([
        "micromamba", "run", "-n", "krepp-tutorial", "gappa", "examine", "heat-tree",
        "--jplace-path", str(JPLACE_OUT), "--mass-norm", "relative",
        "--write-svg-tree", "--svg-tree-shape", "circular",
        "--svg-tree-type", "phylogram", "--svg-tree-ladderize",
        "--svg-tree-stroke-width", "20",
        "--over-color", "#800000",
        "--under-color", "#c8c8cc",
        "--color-list", "#c8c8cc,#d4a0a0,#e07060,#e04040,#d02020,#c01010,#a00808,#800000",
        "--allow-file-overwriting",
        "--out-dir", str(MATERIALS / "results"),
    ], check=True)
    src = MATERIALS / "results" / "tree.svg"
    if not src.exists():
        return
    import shutil, re
    shutil.copy(src, OUT_HEAT_SVG)

    # Thicken high-placement branches based on stroke colour
    svg = OUT_HEAT_SVG.read_text()
    def thicken(m):
        prefix, color, middle = m.group(1), m.group(2), m.group(3)
        r, g, b = int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)
        if r > g * 1.3 and r > b * 1.3:
            w = "160" if r > g * 3 else "80"
            return f'{prefix}stroke="{color}"{middle}stroke-width="{w}"'
        return m.group(0)
    svg = re.sub(r'(<(?:line|path)\b[^>]*?)stroke="(#[0-9a-fA-F]+)"([^/]*?)stroke-width="20"', thicken, svg)
    OUT_HEAT_SVG.write_text(svg)
    print(f"Wrote {OUT_HEAT_SVG}")

    subprocess.run(["rsvg-convert", "-o", str(OUT_HEAT_PNG), "-d", "150", "-w", "3000",
                     str(OUT_HEAT_SVG)], check=True)
    print(f"Wrote {OUT_HEAT_PNG}")


def main():
    run_krepp()
    data = load_mapping()
    for org, mr, ref, d in data:
        print(f"  {org:40s}  mr={mr:5.1f}%  best_ref={ref:25s}  avg_dist={d:.4f}")
    plot_summary(data)
    plot_heat_tree()


if __name__ == "__main__":
    main()

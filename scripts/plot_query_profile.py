#!/usr/bin/env python3
"""
Plot Kraken-style profile.tsv (@@ header): species-level PERCENTAGE bar chart.

Usage:
  python3 scripts/plot_query_profile.py [--profile data-new/profile.tsv] [--out figures/profile_mixture_species.png]
"""

from __future__ import annotations

import argparse
from pathlib import Path


def load_species_rows(path: Path) -> list[tuple[str, float, str]]:
    rows: list[tuple[str, float, str]] = []
    in_table = False
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("@@TAXID"):
                in_table = True
                continue
            if not in_table:
                continue
            if not line.strip():
                break
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            taxid, rank, _taxpath, taxpathsn, pct_s = parts[:5]
            if rank != "species":
                continue
            try:
                pct = float(pct_s)
            except ValueError:
                continue
            label = taxpathsn.split("|")[-1].strip()
            kingdom = taxpathsn.split("|")[0].strip() if "|" in taxpathsn else ""
            rows.append((label, pct, kingdom))
    rows.sort(key=lambda x: x[1])
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", type=Path, default=Path("data-new/profile.tsv"))
    ap.add_argument("--out", type=Path, default=Path("figures/profile_mixture_species.png"))
    args = ap.parse_args()

    species = load_species_rows(args.profile)
    if not species:
        raise SystemExit("No species rows found (expect @@ lines after @@TAXID header).")

    import matplotlib.pyplot as plt

    labels = [s[0] for s in species]
    values = [s[1] for s in species]
    kingdoms = [s[2] for s in species]
    colors = ["#2c5282" if k == "Bacteria" else "#944454" for k in kingdoms]

    fig_h = max(6.0, 0.32 * len(labels) + 1.5)
    fig, ax = plt.subplots(figsize=(10, fig_h), dpi=150)
    y = range(len(labels))
    ax.barh(y, values, color=colors, height=0.72)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Percentage (profile mixture)")
    ax.set_title("Simulated query community: species abundances (profile.tsv)")
    ax.grid(axis="x", alpha=0.35, linestyle="-", linewidth=0.5)
    ax.set_xlim(0, max(values) * 1.08)

    from matplotlib.patches import Patch

    leg = [
        Patch(facecolor="#2c5282", label="Bacteria"),
        Patch(facecolor="#944454", label="Archaea"),
    ]
    ax.legend(handles=leg, loc="lower right", framealpha=0.92)

    plt.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out, bbox_inches="tight")
    plt.close()
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

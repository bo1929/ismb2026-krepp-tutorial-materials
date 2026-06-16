#!/usr/bin/env python3
"""Generate the placement heat-tree figure via gappa.

Requires gappa ≥ 0.9.0 in the krepp-tutorial environment.

Usage:
  python3 scripts/plot_placement_heat_tree.py
"""

from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT.parent / "ismb2026-krepp-tutorial-materials"
JPLACE = MATERIALS / "results" / "placements_default.jplace"
OUT_SVG = ROOT / "figures" / "placement_heat_tree.svg"
OUT_PNG = ROOT / "figures" / "placement_heat_tree.png"


def main():
    if not JPLACE.exists():
        raise FileNotFoundError(
            f"{JPLACE} not found. Run krepp place first."
        )

    # Generate SVG via gappa
    subprocess.run(
        [
            "micromamba", "run", "-n", "krepp-tutorial",
            "gappa", "examine", "heat-tree",
            "--jplace-path", str(JPLACE),
            "--mass-norm", "relative",
            "--write-svg-tree",
            "--svg-tree-shape", "rectangular",
            "--svg-tree-type", "phylogram",
            "--svg-tree-ladderize",
            "--svg-tree-stroke-width", "8",
            "--file-prefix", "heat",
            "--allow-file-overwriting",
            "--out-dir", str(MATERIALS / "results"),
        ],
        check=True,
    )

    # Rename gappa output to heat_tree.svg, then copy to figures
    gappa_out = MATERIALS / "results" / "heattree.svg"
    renamed = MATERIALS / "results" / "heat_tree.svg"
    if gappa_out.exists():
        import shutil
        shutil.move(str(gappa_out), str(renamed))
        shutil.copy(renamed, OUT_SVG)
        print(f"Wrote {OUT_SVG}")

    # Convert to PNG via rsvg-convert
    subprocess.run(
        ["rsvg-convert", "-o", str(OUT_PNG), "-d", "200", str(OUT_SVG)],
        check=True,
    )
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()

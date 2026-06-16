#!/usr/bin/env python3
"""
build_tree.py - build a Newick tree from Mash distances over reference genomes.

Pipeline:
  1. mash sketch (-s 10000 -k 21) from FASTAs listed in input_map.tsv
  2. mash dist (all-vs-all)
  3. neighbor-joining (dendropy) + midpoint rooting

Tip labels match column 1 of input_map.tsv (required for krepp index -t / place).

Usage (from repo root):
  python3 scripts/build_tree.py
  python3 scripts/build_tree.py --input-map data/input_map.tsv -o data/reference_tree.nwk
"""

from __future__ import annotations

import argparse
import io
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional


def resolve_fasta(repo_root: Path, input_map: Path, relpath: str) -> Optional[Path]:
    for base in (input_map.parent, repo_root, repo_root / "data"):
        candidate = (base / relpath).resolve()
        if candidate.exists():
            return candidate
    return None


def load_input_map(path: Path, repo_root: Path) -> list[tuple[str, Path]]:
    """Return (genome_id, absolute_fasta_path) rows in file order."""
    rows: list[tuple[str, Path]] = []
    input_map = path.resolve()
    with open(input_map, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gid, relpath = parts[0].strip(), parts[1].strip()
            fa = resolve_fasta(repo_root.resolve(), input_map, relpath)
            if fa is None:
                sys.exit(f"FASTA missing for {gid}: {relpath} (tried under {input_map.parent})")
            rows.append((gid, fa))
    if not rows:
        sys.exit(f"No genomes in {path}")
    return rows


def mash_sketch_and_dist(fastas: list[Path], mash_bin: str) -> list[tuple[str, str, float]]:
    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        sketch_base = td_path / "refs"
        subprocess.run(
            [mash_bin, "sketch", "-s", "10000", "-k", "21", "-o", str(sketch_base)]
            + [str(p) for p in fastas],
            check=True,
            capture_output=True,
        )
        sketch_path = Path(str(sketch_base) + ".msh")
        cp = subprocess.run(
            [mash_bin, "dist", str(sketch_path), str(sketch_path)],
            check=True,
            capture_output=True,
            text=True,
        )
    rows: list[tuple[str, str, float]] = []
    for line in cp.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        a, b, d = parts[0], parts[1], float(parts[2])
        rows.append((Path(a).name, Path(b).name, d))
    return rows


def stem_to_genome_id(fastas: list[tuple[str, Path]]) -> dict[str, str]:
    """Map Mash sketch labels (file stem) to input_map genome_id."""
    out: dict[str, str] = {}
    for gid, fa in fastas:
        stem = fa.name
        if stem.endswith(".gz"):
            stem = Path(stem[:-3]).stem  # .fna.gz -> accession
        else:
            stem = fa.stem
        out[stem] = gid
        # Mash may report path without .gz in label
        out[fa.stem] = gid
    return out


def distance_matrix(
    genome_ids: list[str],
    mash_rows: list[tuple[str, str, float]],
    stem_map: dict[str, str],
) -> list[list[float]]:
    n = len(genome_ids)
    idx = {gid: i for i, gid in enumerate(genome_ids)}
    dmat = [[0.0] * n for _ in range(n)]
    for a, b, dist in mash_rows:
        ga = stem_map.get(Path(a).stem) or stem_map.get(a)
        gb = stem_map.get(Path(b).stem) or stem_map.get(b)
        if ga is None or gb is None:
            continue
        if ga in idx and gb in idx:
            dmat[idx[ga]][idx[gb]] = dist
    return dmat


def nj_tree_newick(labels: list[str], matrix: list[list[float]]) -> str:
    import dendropy

    buf = io.StringIO()
    buf.write("," + ",".join(labels) + "\n")
    for i, lab in enumerate(labels):
        buf.write(lab + "," + ",".join(f"{matrix[i][j]:.6f}" for j in range(len(labels))) + "\n")
    buf.seek(0)
    pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
        src=buf,
        is_first_row_column_names=True,
        is_first_column_row_names=True,
        delimiter=",",
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
    # Enforce a visible minimum so no branch collapses to zero width in a
    # rectangular phylogram: leaves must terminate horizontally.
    max_depth = 0.0
    for leaf in tree.leaf_node_iter():
        d = 0.0
        n = leaf
        while n.parent_node is not None:
            el = n.edge.length
            d += float(el) if el is not None else 0.0
            n = n.parent_node
        if d > max_depth:
            max_depth = d
    min_edge = max(max_depth * 0.003, 0.0005)
    for node in tree:
        if node.edge.length is not None and node.edge.length < min_edge:
            node.edge.length = min_edge
    out = io.StringIO()
    tree.write(
        file=out,
        schema="newick",
        suppress_rooting=True,
        unquoted_underscores=True,
        suppress_edge_lengths=False,
    )
    return out.getvalue().strip()


def main() -> None:
    ap = argparse.ArgumentParser(description="Build Mash NJ tree from input_map.tsv")
    ap.add_argument(
        "--input-map",
        type=Path,
        default=Path("data/input_map.tsv"),
        help="Two-column TSV: genome_id<TAB>path/to.fasta[.gz]",
    )
    ap.add_argument(
        "-o",
        "--out-tree",
        type=Path,
        default=Path("data/reference_tree.nwk"),
        help="Output Newick path",
    )
    ap.add_argument("--root", type=Path, default=Path("."), help="Repo root for relative FASTA paths")
    ap.add_argument("--mash-bin", default="mash")
    args = ap.parse_args()

    repo = args.root.resolve()
    rows = load_input_map(args.input_map.resolve(), repo)
    genome_ids = [gid for gid, _ in rows]
    fastas = [fa for _, fa in rows]

    sys.stderr.write(f"Mash sketch + dist over {len(genome_ids)} genomes...\n")
    mash_rows = mash_sketch_and_dist(fastas, args.mash_bin)
    stem_map = stem_to_genome_id(rows)
    dmat = distance_matrix(genome_ids, mash_rows, stem_map)

    missing = sum(
        1
        for i in range(len(genome_ids))
        for j in range(i + 1, len(genome_ids))
        if dmat[i][j] == 0.0
    )
    if missing:
        sys.stderr.write(f"warning: {missing} off-diagonal distances are still zero (check label mapping)\n")

    newick = nj_tree_newick(genome_ids, dmat)
    args.out_tree.parent.mkdir(parents=True, exist_ok=True)
    args.out_tree.write_text(newick + "\n", encoding="utf-8")
    sys.stderr.write(f"Wrote {args.out_tree} ({len(genome_ids)} tips)\n")


if __name__ == "__main__":
    main()

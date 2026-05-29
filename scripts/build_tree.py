#!/usr/bin/env python3
"""
build_tree.py - build a Newick tree from mash distances.

Pipeline:
  1. mash sketch each genome (-s 10000 -k 21) into a single .msh
  2. mash dist (all-vs-all) -> square distance matrix
  3. neighbor-joining via dendropy
  4. midpoint root -> tree.nwk

We only use this tree for indexing and placement scaffolding. Branch lengths
encode mash distance approximations; topology is correct at the within-genus
and between-class levels for our marine reference set.
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import dendropy


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes-dir", required=True,
                    help="Directory of <id>.fna files (each id = tip label).")
    ap.add_argument("--keep-ids", required=True,
                    help="File with one tip label per line; tips kept in tree.")
    ap.add_argument("--out-tree", required=True)
    ap.add_argument("--mash-bin", default="mash")
    args = ap.parse_args()

    keep = []
    with open(args.keep_ids) as fh:
        for line in fh:
            line = line.strip()
            if line:
                keep.append(line)
    sys.stderr.write("Building tree over {} tips...\n".format(len(keep)))

    gdir = Path(args.genomes_dir)
    fastas = []
    for gid in keep:
        fa = gdir / f"{gid}.fna"
        if not fa.exists():
            sys.exit(f"FASTA missing: {fa}")
        fastas.append(fa)

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        sketch_path = td / "all.msh"
        sys.stderr.write("  mash sketch...\n")
        run([args.mash_bin, "sketch", "-s", "10000", "-k", "21",
             "-o", str(sketch_path).removesuffix(".msh"), *map(str, fastas)],
            stderr=subprocess.DEVNULL)

        sys.stderr.write("  mash dist (all-vs-all)...\n")
        cp = subprocess.run([args.mash_bin, "dist", str(sketch_path), str(sketch_path)],
                            check=True, capture_output=True, text=True)
        rows = [l for l in cp.stdout.splitlines() if l.strip()]

        def label_of(path):
            return Path(path).stem
        n = len(keep)
        idx = {gid: i for i, gid in enumerate(keep)}
        D = [[0.0] * n for _ in range(n)]
        for line in rows:
            parts = line.split("\t")
            a = label_of(parts[0])
            b = label_of(parts[1])
            d = float(parts[2])
            if a in idx and b in idx:
                D[idx[a]][idx[b]] = d

    matrix_text = "    {}\n".format(n)
    for i, gid in enumerate(keep):
        name = gid[:10].ljust(10)
        matrix_text += name + " " + " ".join("{:.6f}".format(D[i][j]) for j in range(n)) + "\n"

    pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
        src=_to_csv(keep, D), is_first_row_column_names=True, is_first_column_row_names=True,
        delimiter=",")
    nj = pdm.nj_tree()
    nj.encode_bipartitions()
    nj.reroot_at_midpoint(update_bipartitions=False)

    nj.write(path=args.out_tree, schema="newick", suppress_rooting=True,
             unquoted_underscores=True, suppress_edge_lengths=False)
    sys.stderr.write("Wrote {}\n".format(args.out_tree))


def _to_csv(labels, matrix):
    import io
    buf = io.StringIO()
    buf.write("," + ",".join(labels) + "\n")
    for i, lab in enumerate(labels):
        buf.write(lab + "," + ",".join("{:.6f}".format(matrix[i][j]) for j in range(len(labels))) + "\n")
    buf.seek(0)
    return buf


if __name__ == "__main__":
    main()

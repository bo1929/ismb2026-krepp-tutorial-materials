#!/usr/bin/env python3
"""
Simulate a mixed-community Illumina library with ART, proportional to species
abundances in profile.tsv (species rows). Concatenates ART outputs into one FASTQ.

Single-end (ART):  art_illumina -ss HS25 -l L -f c -na -i INPUT.fna -o PREFIX
For paired-end with insert stdev -s, ART also requires -p -m (mean fragment).

Default: 1,000,000 single-end reads, read length 150 (HiSeq 2500).
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def parse_species_profile(path: Path) -> list[tuple[str, float]]:
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("@") or line.startswith("#"):
                continue
            if line.startswith("@@"):
                continue
            parts = line.split("\t")
            if len(parts) < 5 or parts[1] != "species":
                continue
            taxid = parts[0].strip()
            pct = float(parts[4])
            rows.append((taxid, pct))
    return rows


def fasta_total_length(fa: Path) -> int:
    n = 0
    with open(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            n += len(line.strip())
    return n


def allocate_reads(weights: list[float], total: int) -> list[int]:
    s = sum(weights)
    if s <= 0:
        raise SystemExit("abundance sum must be positive")
    raw = [total * w / s for w in weights]
    out = [int(x) for x in raw]
    rem = total - sum(out)
    if rem <= 0:
        return out
    frac = sorted(range(len(weights)), key=lambda i: raw[i] - out[i], reverse=True)
    for i in frac[:rem]:
        out[i] += 1
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", type=Path, default=Path("profile.tsv"))
    ap.add_argument("--genomes-dir", type=Path, default=Path("profile_genomes"))
    ap.add_argument("-o", "--output", type=Path, default=Path("query_mixture.fq"))
    ap.add_argument("-n", "--total-reads", type=int, default=1_000_000)
    ap.add_argument("-l", "--read-length", type=int, default=150,
                    help="ART -l (HiSeq 2500: 125 or 150)")
    ap.add_argument("--art", default="art_illumina", help="ART binary name or path")
    args = ap.parse_args()

    if shutil.which(args.art) is None:
        sys.stderr.write(f"not found on PATH: {args.art}\n")
        sys.exit(1)

    spec = parse_species_profile(args.profile)
    if not spec:
        sys.stderr.write(f"no species rows in {args.profile}\n")
        sys.exit(1)

    taxids = [t for t, _ in spec]
    weights = [w for _, w in spec]
    nreads = allocate_reads(weights, args.total_reads)
    if sum(nreads) != args.total_reads:
        sys.stderr.write("internal error: read allocation does not sum to target\n")
        sys.exit(1)

    tmp = Path(tempfile.mkdtemp(prefix="art_mix_"))
    try:
        pieces: list[Path] = []
        for (taxid, pct), nr in zip(spec, nreads):
            if nr <= 0:
                continue
            fa = args.genomes_dir / f"{taxid}.fna"
            if not fa.is_file():
                sys.stderr.write(f"missing genome: {fa}\n")
                sys.exit(1)
            glen = fasta_total_length(fa)
            if glen <= 0:
                sys.stderr.write(f"empty fasta: {fa}\n")
                sys.exit(1)
            fold = (nr * args.read_length) / glen
            prefix = tmp / taxid
            cmd = [
                args.art,
                "-ss", "HS25",
                "-l", str(args.read_length),
                "-f", str(fold),
                "-na",
                "-i", str(fa),
                "-o", str(prefix),
            ]
            sys.stderr.write(
                f"{taxid}  abundance={pct:.4g}%  reads={nr}  len={glen}  -f={fold:.6g}\n"
            )
            cp = subprocess.run(cmd, capture_output=True, text=True)
            if cp.returncode != 0:
                sys.stderr.write(cp.stderr or cp.stdout or "ART failed\n")
                sys.exit(cp.returncode)
            fq = Path(str(prefix) + ".fq")
            if not fq.is_file():
                sys.stderr.write(f"expected ART output: {fq}\n")
                sys.exit(1)
            pieces.append(fq)

        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "wb") as out:
            for fq in pieces:
                with open(fq, "rb") as inh:
                    shutil.copyfileobj(inh, out)

        # ART -f yields slightly fewer reads than (fold*genome_len/read_length); top up.
        with open(args.output, "rb") as fh:
            n_lines = sum(1 for _ in fh)
        have = n_lines // 4
        need = args.total_reads - have
        if need > 0:
            top_taxid = max(zip(taxids, nreads), key=lambda x: x[1])[0]
            fa = args.genomes_dir / f"{top_taxid}.fna"
            fold = (need * args.read_length) / fasta_total_length(fa) * 2.0
            corr = tmp / "topup"
            cmd = [
                args.art,
                "-ss", "HS25",
                "-l", str(args.read_length),
                "-f", str(fold),
                "-na",
                "-i", str(fa),
                "-o", str(corr),
            ]
            cp = subprocess.run(cmd, capture_output=True, text=True)
            if cp.returncode != 0:
                sys.stderr.write(cp.stderr or "ART top-up failed\n")
                sys.exit(cp.returncode)
            fq_top = Path(str(corr) + ".fq")
            n_need_lines = need * 4
            with open(args.output, "ab") as out, open(fq_top, "rb") as inh:
                for i in range(n_need_lines):
                    line = inh.readline()
                    if not line:
                        sys.stderr.write(
                            f"warning: top-up FASTQ ended after {i // 4} reads\n"
                        )
                        break
                    out.write(line)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    n_lines = sum(1 for _ in open(args.output, "rb"))
    sys.stderr.write(
        f"Wrote {args.output}  ({n_lines} lines, "
        f"{n_lines // 4} reads)\n"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
simulate_queries.py - generate FASTQ query reads with a realistic abundance profile.

Pure-Python read simulator (no wgsim dependency). Produces single-end 150bp reads
with a configurable per-base substitution error rate. Deterministic given a seed.

Usage:
    python3 simulate_queries.py \
        --profile abundance_profile.tsv \
        --genomes-dir data/genomes \
        --metadata data/metadata.tsv \
        --total-reads 5000 \
        --out data/query.fq

Each profile line is: <short_id> <TAB> <relative abundance>. Counts are floor()
of total * abundance with leftover assigned to the most abundant source so the
sum matches --total-reads. Read IDs are formatted as
    @<role>:<short_id>:<idx>
where role is 'novel' for holdout genomes and 'tip' for in-reference ones.
"""

import argparse
import gzip
import os
import random
import sys
from pathlib import Path

ALPHABET = "ACGT"
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    open_fn = gzip.open if str(path).endswith(".gz") else open
    seqs = []
    header = None
    chunks = []
    with open_fn(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(chunks).upper()))
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            seqs.append((header, "".join(chunks).upper()))
    return seqs


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def mutate(seq, rng, error_rate):
    if error_rate <= 0:
        return seq
    out = list(seq)
    for i, b in enumerate(out):
        if b not in ALPHABET:
            continue
        if rng.random() < error_rate:
            choices = [x for x in ALPHABET if x != b]
            out[i] = rng.choice(choices)
    return "".join(out)


def sample_read(contigs, rng, read_len, error_rate):
    total = sum(len(s) for _, s in contigs)
    if total < read_len:
        return None
    for _ in range(50):
        target = rng.randrange(total - read_len)
        offset = 0
        for _, s in contigs:
            if target < offset + len(s) - read_len:
                local = target - offset
                window = s[local:local + read_len]
                if "N" in window:
                    break
                if rng.random() < 0.5:
                    window = revcomp(window)
                return mutate(window, rng, error_rate)
            offset += len(s)
    return None


def load_profile(path):
    items = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            items.append((parts[0], float(parts[1])))
    total = sum(a for _, a in items)
    if total <= 0:
        sys.exit("Empty abundance profile.")
    return [(g, a / total) for g, a in items]


def load_roles(metadata_path):
    roles = {}
    with open(metadata_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            row = dict(zip(header, parts))
            roles[row["genome"]] = row.get("role", "reference")
    return roles


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True)
    ap.add_argument("--genomes-dir", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--total-reads", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--error-rate", type=float, default=0.005)
    args = ap.parse_args()

    profile = load_profile(args.profile)
    roles = load_roles(args.metadata)
    rng = random.Random(args.seed)

    counts = [(gid, max(1, int(args.total_reads * w))) for gid, w in profile]
    deficit = args.total_reads - sum(c for _, c in counts)
    counts[0] = (counts[0][0], counts[0][1] + deficit)

    sys.stderr.write("Simulating {} reads from {} sources:\n".format(args.total_reads, len(counts)))
    written = 0
    truth_path = Path(args.out).with_suffix(".truth.tsv")
    with open(args.out, "w") as out, open(truth_path, "w") as truth:
        truth.write("read_id\trole\tsource_genome\tabundance\tcount\n")
        for gid, n in counts:
            role = "novel" if roles.get(gid) == "holdout" else "tip"
            fa = Path(args.genomes_dir) / f"{gid}.fna"
            contigs = read_fasta(fa)
            produced = 0
            attempts = 0
            ab = next(w for g, w in profile if g == gid)
            while produced < n and attempts < n * 10:
                read = sample_read(contigs, rng, args.read_len, args.error_rate)
                attempts += 1
                if read is None:
                    continue
                rid = f"{role}_{gid}_{produced}"
                out.write(f"@{rid}\n{read}\n+\n{'I' * len(read)}\n")
                truth.write(f"{rid}\t{role}\t{gid}\t{ab:.4f}\t{n}\n")
                produced += 1
                written += 1
            sys.stderr.write(f"  {role:<8s} {gid:<22s} {produced:>5d} reads ({ab*100:5.2f}%)\n")
    sys.stderr.write("Wrote {} reads to {} (ground truth: {})\n".format(written, args.out, truth_path))


if __name__ == "__main__":
    main()

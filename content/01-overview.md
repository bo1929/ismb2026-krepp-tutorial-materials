# Overview

**krepp** is a *k*-mer-based tool for searching and finding query sequences in large reference genome collections, estimating the distance between reads and reference genomes, and placing them on backbone phylogenies and/or taxonomies.
It uses locality-sensitive hashing (LSH), a phylogeny-guided *k*-mer coloring algorithm, and a maximum-likelihood framework to achieve these.

---

## What we will cover

This tutorial walks through the core krepp workflow using a controlled mock microbial community consisting of 20 marine genomes simulated at 100,000 Illumina reads.
We will work with two small reference sets to avoid heavy computation during the tutorial.
We will index one toy reference set with controlled novelty and limited relatedness to queries and directly download another, diverse but heavily down-sampled microbial index.

The tutorial consists of five main sections:

1. Installation of the tools and downloading the tutorial datasets.
2. Indexing reference genomes and using one of the available pre-built indexes.
3. Estimating read-to-genome distances and interpretation of distances.
4. Placing reads on a reference phylogeny (or a taxonomy) and visualising the result with `gappa`.
5. Exercises with a more realistic reference index and conclusions.

## Resources

- [Slides](https://bo1929.github.io/presentations/krepp-ismb-2026.pdf) for the first part of this tutorial.
- [The main paper (Genome Biology)](https://doi.org/10.1186/s13059-026-03999-y) with methodological details and extensive benchmarking.
- [Another tutorial (in IMSI, 2025)](https://bo1929.github.io/presentations/krepp-imsi-2025.pdf) for an earlier version which remains mostly valid.
- [Slides (from RECOMB 2025)](https://bo1929.github.io/presentations/krepp-recomb-2025.pdf) focusing on method and results.

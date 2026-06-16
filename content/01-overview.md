# Overview

## What is krepp?

**krepp** is a *k*-mer-based tool for searching and finding query sequences in large reference genome collections, estimating the distance between them, and placement of on backbone phylogenies and/or taxonomies. It uses locality-sensitive hashing (LSH), a phylogeny-guided *k*-mer coloring algorithm, and a maximum-likelihood framework to achieve these.

---

## What we will cover

This tutorial walks through the core krepp workflow using a controlled mock community of 20 marine microbial genomes simulated at 100,000 Illumina reads. We will work with two reference sets:

- **A 31-genome toy index** with controlled taxonomic novelty, for learning and interpretation.
- **A Web of Life (WoL) tiny index** (>10,000 genomes), to demonstrate the same pipeline at scale.

| Step | Command | What it does |
|------|---------|---------------|
| 1. Index | `krepp index` | Build a searchable index from reference genomes |
| 2. Distance | `krepp dist` | Estimate per-read distances to references |
| 3. Placement | `krepp place -t` | Place reads on a reference phylogeny |
| 4. Taxonomy | `krepp place -l` | Place reads on a taxonomic tree via lineages |
| 5. Visualize | gappa | Render placements as heat-trees |

1. Build a krepp index from reference genomes
2. Estimate read-to-reference distances and interpret them by novelty level.
3. Place reads on a reference phylogeny and visualise the result with gappa.
4. Use taxonomic lineages to annotate placements and compare against ground truth.
5. Apply the same commands to a large-scale index (WoL) and compare results.

## Resources

- Slides for the first part of this tutorial (key results and overview)
- Genome Biology paper (methodological details and detailed benchmarking)
- IMSI Software Workshop (an earlier tutorial version but remains valid)
- RECOMB 2025 (a presentation focusing on method and results)

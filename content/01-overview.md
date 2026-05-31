# Overview

## Background

## What krepp can do

| Command | Input | Output |
|---------|-------|--------|
| `krepp index` | reference FASTAs + tree | color-index |
| `krepp inspect` | index | index statistics |
| `krepp dist` | index + reads | per-read k-mer distances |
| `krepp dist --summarize` | index + reads | per-reference abundance table |
| `krepp place` | index + reads | phylogenetic placements (jplace) |
| `krepp place --tabular` | index + reads | tabular placements (TSV) |
| `krepp place --summarize` | index + reads | per-branch abundance table |
| `krepp place -l lineages.tsv` | index + reads + lineages | taxonomic placements |


## Learning outcomes

1. Build a krepp index from reference genomes and a phylogenetic tree.
2. Use `krepp dist` to estimate read-to-reference distances and infer
   community composition.
3. Use `krepp place` to place reads on a tree and identify the branch
   accumulating novel-organism reads.
4. Add taxonomic lineages to placement output and interpret the results.
5. Visualize placements as a heat tree using gappa.

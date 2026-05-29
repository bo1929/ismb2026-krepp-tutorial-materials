# Overview

## Background

Shotgun metagenomic reads are short (100-300 bp) and often diverge from any
genome in your database. Classic alignment-based methods require a read to map
cleanly to at least one reference. For highly divergent or truly novel
organisms, alignment fails silently: the read is simply discarded.

**krepp** takes a similar route. It works with *k-mers* - exact subsequences
of length k - and computes an evolutionary distance between a read and every
reference without requiring a global alignment.

---

## Core ideas

### *k*-mer distance

Two sequences are similar if they share many k-mers. A Jukes-Cantor correction
converts the fraction of shared k-mers into an evolutionary distance:

```
d_JC = -3/4 * ln(1 - 4/3 * (1 - jaccard(k-mers_read, k-mers_ref)))
```

For a 150 bp read with k = 27, there are up to 124 minimizer k-mers. Even at
5% divergence (typical within-genus) many k-mers are conserved, giving a
reliable signal.

### Locality-sensitive hashing (LSH)

Storing all k-mers of all references verbatim would require tens of GB. krepp
compresses the k-mer space with **locality-sensitive hashing**: k-mers are
projected into a smaller space such that similar k-mers hash to the same
bucket with high probability. Combined with **minimizers** (a sliding window
retaining only the lexicographically smallest k-mer), the index shrinks
dramatically while retaining statistical power.

### Maximum pseudo-likelihood placement

Once per-reference distances are computed for a read, krepp evaluates *where*
on the reference tree the read most plausibly evolved. It assigns
**likelihood-weighted ratios (LWR)** to every branch, outputting a jplace
file or a tabular placement table.

---

## What krepp can do

| Command | Input | Output |
|---------|-------|--------|
| `krepp index` | reference FASTAs + tree | color-index |
| `krepp inspect` | color-index | index statistics |
| `krepp dist` | index + reads | per-read k-mer distances |
| `krepp dist --summarize` | index + reads | per-reference abundance table |
| `krepp place` | index + reads | phylogenetic placements (jplace) |
| `krepp place --tabular` | index + reads | tabular placements (TSV) |
| `krepp place --summarize` | index + reads | per-branch abundance table |
| `krepp place -l lineages.tsv` | index + reads + lineages | taxonomic placements |

---

## Our biological scenario

!!! abstract "The Challenge"
    You have **100 000** shotgun reads (`data/query_mixture.fq.gz`) drawn
    from a synthetic mixture of **20 profile species** (marine and related taxa).
    A compact reference genome set of **31 assemblies** is indexed in krepp (see
    `data/input_map.tsv`). **`data/reference_info.tsv`** describes the designed
    relationship between the query organisms and the reference genome set.
    Which lineages explain the reads? Where do reads land when the closest reference is
    another species or genus?

The references mix **same-species** alternatives, **same-genus** surrogates,
surrogates with coarse shared ranks (e.g. `kingdom`, `family`), panel extras, and **class-disjoint** references,
documented in `data/reference_info.tsv`. Several mixture
taxa therefore hit references that are phylogenetically close but **not** the
same genome as in `data/query_info.tsv` - the textbook signature for
**nearest-neighbor** abundance inflation and **internal-branch** placements.

---

## Learning outcomes

1. Build a krepp color-index from reference genomes and a phylogenetic tree.
2. Use `krepp dist` to estimate read-to-reference distances and infer
   community composition.
3. Use `krepp place` to place reads on a tree and identify the branch
   accumulating novel-organism reads.
4. Add taxonomic lineages to placement output and interpret the results.
5. Visualize placements as a heat tree using gappa.

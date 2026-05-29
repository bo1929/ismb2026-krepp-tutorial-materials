# Synthesis

## What did we discover?

| Step | Tool | Key finding |
|------|------|-------------|
| Community composition | `krepp dist --summarize` | References phylogenetically closest to mixture genomes absorb extra mass when assemblies differ |
| Per-read distances | `krepp dist` | Non-zero distances to the best neighbor signal strain or genus mismatch |
| Phylogenetic placement | `krepp place --tabular` | Terminal vs internal (`DISTAL_NODE = NA`) splits quantify unresolved divergence |
| Taxonomic placement | `krepp place -l lineages.tsv` | Rank labels contextualize those splits once lineages are supplied |
| Visualization | gappa heat tree | Thick surrogate branches mirror summarize inflation |

!!! success "Conclusion"
    The `data-new` mixture pairs **20 profile species** with a deliberately
    structured **31-genome** index (`input_map.tsv`). **`reference_info.tsv`**
    documents **35** panel rows (five assemblies are described there but omitted
    from the index - see Dataset). Expect abundance summaries and
    placements to emphasize **nearest indexed relatives**, not necessarily the
    query assembly accession listed in `query_info.tsv`.

---

## The novel organism workflow

```
1. dist --summarize  ->  flag references with unexpectedly high mass vs profile
2. place --tabular   ->  look for NA distal nodes (internal edges)
3. place -l lineages ->  attach rank names once lineages exist
4. gappa heat-tree   ->  visualize concentrated branches
```

---

## krepp vs. alignment-based tools

| Property | krepp | bowtie2 + MetaPhlAn |
|----------|-------|---------------------|
| Detects novel taxa | Yes (internal branch signal) | No (unmapped read) |
| Speed | Fast (k-mer lookup) | Moderate (alignment) |
| Novel classification rank | Any level | None |
| Reference requirement | Phylogenetic tree | Database only |
| Output format | jplace / TSV | TSV |

---

## Exercises

??? question "Exercise 1 - Abundance reconciliation"
    Compare `krepp dist --summarize` with species percentages in `profile.tsv`.
    Which genome ids reconcile only after summing multiple strains or surrogates?

??? question "Exercise 2 - Panel edits"
    Remove one reference row from `input_map.tsv` (and rebuild the index with a
    matching pruned tree). How do placements for reads tied to that taxon change?

??? question "Exercise 3 - Distance threshold sensitivity"
    Run `krepp dist --filter --no-multi` and count how many reads still report a
    hit for a `novel` mixture genome. Repeat with `--tau 1` and `--tau 10`.

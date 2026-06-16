# Conclusions

We covered the core krepp workflow: indexing, distance estimation, phylogenetic placement, and taxonomic placement using a small toy reference set.

**Distance estimation** (`krepp dist`) computes per-read Hamming distances to reference genomes, akin to the ratio of mismatches one would obtain through alignment, but sensitive to more divergent sequences.
The default command produces a tab-separated table of read-to-reference mappings and maximum likelihood distances.
Filtering (`--filter`) retains only statistically close hits; the `--summarize` flag aggregates results into operational genomic unit counts, which can be used for downstream analysis and sample comparison.

**Phylogenetic placement** (`krepp place -t`) finds placements of query reads on a given reference phylogeny, producing a `.jplace` file by default, a format that can be further processed by tools such as gappa, iTOL, and QIIME 2.
**Taxonomic placement** (`krepp place -l`) uses a lineage file to assign reads to named clades.
Both support tabular (`--tabular`) and summarised (`--summarize`) outputs, similar to distance estimation.

---

## Exercises with Web of Life

For further exercises, you can use the same type of commands with a slightly larger index we downloaded earlier.
This is the tiny (i.e., heavily minimized and downsampled) version of the Web of Life index (>10,000 microbial genomes, see [this page](https://github.com/bo1929/krepp/wiki/Available-reference-indexes) for the full versions) on the same query mixture:

<img src="figures/wol_query_summary.png" alt="Per-query mapping rate and best-reference distance for WoL index." style="display: block; margin: 1rem auto; width: 100%; max-width: 850px; height: auto;" />

**Figure 6.** Per-query mapping rate and mean distance to the best-matching reference, using the WoL index.
With denser sampling and a more diverse representation, mapping rates are uniformly high and best-reference distances are consistently low, meaning that all these queries have relatively close matches in this index compared to the toy index we built.

<img src="figures/wol_heat_tree-e.png" alt="Heat-tree of placement mass on the WoL phylogeny." style="display: block; margin: 1rem auto; width: 100%; max-width: 900px; height: auto;" />

**Figure 7.** Placement heat-tree on the WoL reference phylogeny (>10,000 leaves). Since the toy sample we generated has low complexity, only a few clades/branches are highlighted. Real and complex samples could result in more widely spread placements.

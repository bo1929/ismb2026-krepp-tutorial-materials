# Dataset

We use one controlled mock metagenome: **100,000 simulated Illumina reads**
from **20 query genomes**. The query profile is known, so every later distance,
placement, and abundance result can be checked against the organisms that truly
generated the reads.

The tutorial uses two reference contexts. First, a **31-genome toy reference
genome set** lets us inspect the behavior by hand. Second, we install the
**WoL-v1 tiny index**, a Web of Life reference index with roughly **10,000
microbial references** and a matching phylogeny, under
`data/index-WoLv1-tiny/`. The toy set is for interpretation; WoL-v1 is for
showing the same workflow on a broader database.

The toy reference genome set is intentionally uneven. Each query and reference
genome has a **`role`** in `query_info.tsv` / `reference_info.tsv`: the
**deepest shared taxonomic rank** (species, genus, family, ...) between that
genome and its closest partner in the other set. Ranks come from
`query_taxonomy.tsv` (Kraken-style paths) and NCBI lineages for references;
recompute with `python3 scripts/update_taxonomy_roles.py`.

Some queries match a reference at **species** rank (same taxon, different
assembly). Others share only **genus**, **family**, or **kingdom** with their
best reference — the distant cases that drive placement on internal branches.

<img src="figures/dataset_profile.png" alt="Query profile composition: kingdom bar, class bar, and per-species abundance chart." style="display: block; margin: 1rem auto; width: 100%; max-width: 760px; height: auto;" />

**Figure 1.** Query profile composition. Top: kingdom- and class-level
stacked bars showing the overall taxonomic spread. Bottom: all 20 species
with their relative abundances, coloured by class.

<img src="figures/dataset_reference.png" alt="How query reads relate to the reference panel, and the reference panel composition by role." style="display: block; margin: 1rem auto; width: 100%; max-width: 720px; height: auto;" />

**Figure 2.** Left: query read abundance grouped by each species' **`role`**
(deepest shared rank with the best-matching reference). Right: the 31-genome
reference panel counted by the same rank labels. Coarser ranks (e.g. genus,
family, kingdom) indicate a more distant best reference.

---

## Files

| File or directory | What it contains |
|-------------------|------------------|
| `data/profile.tsv` | Ground-truth taxonomic profile and abundances across ranks. |
| `data/query_taxonomy.tsv` | Taxonomy paths (`|`-separated) for each query taxon. |
| `data/query_info.tsv` | The 20 query organisms; `role` = match rank to best reference. |
| `data/reference_taxonomy.tsv` | Cached NCBI taxid paths for reference genomes. |
| `data/query_genomes/` | Genome assemblies for those query organisms. |
| `data/query_mixture.fq.gz` | The simulated read mixture used in the tutorial. |
| `data/input_map.tsv` | The 31 reference genomes used to build the small toy index. |
| `data/reference_tree.nwk` | Mash distance NJ tree over the references (tip labels match `input_map` IDs). |
| `data/reference_info.tsv` | Each reference genome; `role` = match rank to best query. |
| `data/index-WoLv1-tiny/` | The installed WoL-v1 tiny index used as a larger reference context. |

---

## What to look for

Good results should recover the abundant close matches, keep same-genus reads
near the right neighborhood, and place the distant fraction on plausible
internal branches instead of treating every read as an exact reference hit.

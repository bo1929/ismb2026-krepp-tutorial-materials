# Taxonomic Placement - `krepp place -l`

Phylogenetic placement tells you *where* on a tree a read belongs.
Taxonomic placement adds the *name*: kingdom, phylum, class, order, family,
genus, or species.

krepp integrates taxonomy via a lineage file that annotates every **indexed**
genome id (column 1 of `data-new/input_map.tsv`) with a GTDB-style path.

!!! warning "Build `lineages_refs.tsv` for data-new"
    The **`data-new/`** bundle does not ship `lineages_refs.tsv`. Export ranks from
    GTDB/NCBI for each genome id, one tab-separated lineage string per row:

    ```
    Alteromonas_macleodii   k__Bacteria; p__...; g__Alteromonas; s__...
    ```

    Names must match the tree leaves exactly.

---

## The lineage file

```bash
head data-new/lineages_refs.tsv
```

```
Alteromonas_macleodii   k__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;...
Pyrodictium_occultum    k__Archaea;p__Crenarchaeota;...
```

Format: `<genome_id><TAB><semicolon-separated lineage>` with rank prefixes
(`k__`, `p__`, `c__`, `o__`, `f__`, `g__`, `s__`).

---

## Running taxonomic placement

```bash
krepp place \
    -i data-new/ref_index \
    -q data-new/query_mixture.fq.gz \
    -l data-new/lineages_refs.tsv \
    --tabular \
    --num-threads 4 \
    > results/place_tax.tsv
```

The header comment embeds the taxonomy-expanded tree. Distal nodes can repeat the
short genome id and longer Latin labels attached to the same edge geometry.

---

## Interpreting novel read placements

```bash
grep '^novel_' results/place_tax.tsv | awk '{print $2}' | sort | uniq -c | sort -rn | head
```

Reads originating from mixture genomes without **exact** assembly matches often
still classify within the expected **genus or family** because surrogates occupy
neighboring leaves.

!!! tip "Rank-level novelty detection"
    If an organism is novel at the genus level, reads tend to pool on internal
    branches shared by references of that genus. Higher-rank novelty moves mass
    toward deeper internal edges in the taxonomy-decorated view.

??? question "Q 6.1"
    If a profile species were missing entirely from the panel but its **family**
    were represented, where would you expect placement density to accumulate?

??? question "Q 6.2"
    When `DISTAL_NODE` lists both a short genome id and a Latin species label on
    adjacent edges with split LWR, what does that say about minimizer distances
    to each decorated node?

---

## Summarize at taxonomic level

```bash
krepp place -i data-new/ref_index -q data-new/query_mixture.fq.gz \
    -l data-new/lineages_refs.tsv \
    --summarize --num-threads 4 > results/place_tax_summarize.tsv
```

This produces per-taxonomic-node abundance estimates across the expanded tree.

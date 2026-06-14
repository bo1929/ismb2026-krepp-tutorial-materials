# Visualization - gappa Heat Trees

The `.jplace` file from `krepp place` can be passed directly to **gappa**
to produce a heat tree where branch width and color reflect placement density.

---

## Install gappa

```bash
conda install -c conda-forge -c bioconda gappa -y
gappa --version
```

---

## Generate a heat tree

```bash
mkdir -p results/heat_tree

gappa examine heat-tree \
    --jplace-path results/placements.jplace \
    --out-dir      results/heat_tree/ \
    --svg-tree-width 1000
```

gappa writes `results/heat_tree/heat_tree.svg`. Open in any browser or SVG
viewer.

---

## Reading the heat tree

- **Thick blue branches** - many placements; dominant organisms in the sample.
- **Thin / grey branches** - few or no placements.

Branches leading to the **closest available references** (different assembly,
species, or genus than the mixture source) often appear thicker than naive
abundance-from-profile estimates would predict - parallel to
`krepp dist --summarize` inflation.

Internal branches above pairs of related references thicken when reads distribute
between terminals and the ancestral edge.

---

## Abundance-weighted mode

```bash
gappa examine heat-tree \
    --jplace-path    results/placements.jplace \
    --out-dir        results/heat_tree/ \
    --mass-policy    accumulated \
    --svg-tree-width 1000
```

With `--mass-policy accumulated`, color represents cumulative LWR mass per
branch rather than raw read count, better reflecting relative allocation across
the tree.

---

## Separate tip vs. novel reads

```bash
# Novel reads only (simulator prefix)
grep '^novel_' results/place_tab.tsv > results/novel_reads_place.tsv

# Tip reads only
grep -v '^novel_\|^SEQ_ID\|^#' results/place_tab.tsv > results/tip_reads_place.tsv
```

Novel-only heat maps highlight internal versus terminal mass for rows flagged in
`reference_info.tsv`.

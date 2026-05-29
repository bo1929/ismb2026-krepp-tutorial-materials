# Phylogenetic Placement - `krepp place`

`krepp place` extends distance estimation to tree placement: given the
per-reference distances for a read, it computes a likelihood-weighted
placement on every branch of the reference phylogeny.

Commands assume **`data-new/ref_index`** and **`data-new/query_mixture.fq.gz`**.

---

## jplace output

```bash
krepp place -i data-new/ref_index -q data-new/query_mixture.fq.gz \
    --num-threads 4 -o results/placements.jplace
```

The `.jplace` format (JSON) is the standard for evolutionary placement
algorithms (EPA, pplacer). It embeds the reference tree and a placement
record per read, and is directly consumable by gappa and iTOL.

---

## Tabular placement output

```bash
krepp place -i data-new/ref_index -q data-new/query_mixture.fq.gz \
    --tabular --num-threads 4 > results/place_tab.tsv
```

Output columns: `SEQ_ID`, `DISTAL_NODE`, `EDGE_NUM`, `LWR`, `DIST`

```
SEQ_ID                   DISTAL_NODE              EDGE_NUM  LWR     DIST
tip_28108_0              Alteromonas_macleodii    12        1.0000  0.0091
novel_54248_0            Pyrodictium_occultum     45        0.6200  0.0382
novel_54248_0            NA                       46        0.3800  0.0599
```

Edge integers depend on your Newick; treat them as opaque IDs stable within one run.

!!! danger "Internal branch = unresolved divergence"
    `DISTAL_NODE = NA` marks mass on an **internal** branch - reads evolve from
    an ancestor **between** named references. That is the hallmark placement when
    the sample contains variation **not** represented on either adjacent leaf
    (different assembly, species surrogate, or broader novelty).

Filter reads that were simulated as holdouts:

```bash
grep '^novel_' results/place_tab.tsv | head -10
```

Compare placement concentration on surrogate terminals vs internal branches for
rows whose **`role`** in `reference_info.tsv` is coarser than `species` (e.g.
`genus`, `family`, `kingdom`).

??? question "Q 5.1"
    Sketch the subtree around `Alteromonas_macleodii` and its genus neighbors in
    the reference tree. Where should reads from the profile genome accumulate if
    the indexed reference is a **different species** in *Alteromonas*?

??? question "Q 5.2"
    What does LWR = 1.0 on a terminal genome id mean? What does a split LWR
    across a terminal edge and an internal edge imply about residual divergence?

---

## Community-level placement - `--summarize`

```bash
krepp place -i data-new/ref_index -q data-new/query_mixture.fq.gz \
    --summarize --num-threads 4 > results/place_summarize.tsv
```

Output columns: `DISTAL_NODE`, `EDGE_NUM`, `WEIGHTED_COUNT`, `SEQUENCE_ABUNDANCE`

Internal-branch rows (`DISTAL_NODE = NA`) aggregate uncertainty between references.
Surrogate references with a coarse shared `role` often deposit measurable mass there even
when `--summarize` also peaks on the closest labeled genome.

??? question "Q 5.3"
    Why can internal edges carry placement mass even when every genome in the
    mixture is represented **somewhere** in the reference genome set?

---

## Sensitivity controls

| Flag | Effect |
|------|--------|
| `--tau 3` | Stricter cutoff - fewer noisy placements |
| `--tau 10` | Looser - useful for very divergent novel taxa |
| `--no-multi` | Single best placement per read only |
| `--no-filter` | Disable distance pre-filter; include all branches |

# Distance Estimation - `krepp dist`

`krepp dist` computes the Jukes-Cantor k-mer distance between every read and
every reference genome in the index. It is the fastest krepp operation and the
most direct way to profile community composition.

Examples below use **`data-new/ref_index`** and **`data-new/query_mixture.fq.gz`**
(or the decompressed `.fq`). Create `results/` as needed.

---

## Per-read distances (default)

```bash
krepp dist -i data-new/ref_index -q data-new/query_mixture.fq.gz --num-threads 4 \
    > results/dist_default.tsv
```

Output columns: `SEQ_ID`, `REFERENCE_NAME`, `DIST`

```
# software: krepp  version: v0.8.0
SEQ_ID                   REFERENCE_NAME           DIST
tip_28108_0              Alteromonas_macleodii    0.0091
novel_54248_0            Pyrodictium_occultum     0.0382
novel_54248_0            Thermus_thermophilus     0.1124
```

Illustrative only - your numeric distances depend on the exact tree and index.

Each read can match multiple references within the reporting threshold. Reads
from a mixture genome whose assembly **differs** from the closest indexed
reference often land at **non-zero** distance on that neighbor - the raw signal
for downstream placement.

---

## Best-hit filter

`--filter` retains only the single closest reference per read:

```bash
krepp dist -i data-new/ref_index -q data-new/query_mixture.fq.gz --filter \
    --num-threads 4 > results/dist_filter.tsv
```

??? question "Q 4.1"
    For a read simulated from a query whose `role` in `query_info.tsv` is coarse
    (e.g. `kingdom` or `family`), which reference label tends to win `--filter`?
    Does the mash / k-mer distance match your expectation from the phylogeny?

---

## Community abundance - `--summarize`

For a quick abundance estimate without storing per-read output:

```bash
krepp dist -i data-new/ref_index -q data-new/query_mixture.fq.gz \
    --summarize --num-threads 4 > results/dist_summarize.tsv
```

Output columns: `REFERENCE_NAME`, `WEIGHTED_COUNT`, `SEQUENCE_ABUNDANCE`

```
REFERENCE_NAME           WEIGHTED_COUNT   SEQUENCE_ABUNDANCE
Spiribacter_roseus       ...
Nitrosococcus_watsonii   ...
Alteromonas_macleodii    ...
```

!!! danger "Nearest-neighbor inflation"
    When mixture reads come from assembly **A** but the index holds a related
    assembly **B** (another strain, species in the same genus, or a deliberate
    surrogate with a coarser shared rank, e.g. `genus` or `kingdom`), **`krepp dist --summarize`**
    credits much of that mass to **B**. Compare summaries to `profile.tsv` and
    `reference_info.tsv` to see which references act as sinks.

??? question "Q 4.2"
    Without opening `profile.tsv`, which references would you inspect first as
    likely sinks for reads from **different assemblies** of the same taxon?

??? question "Q 4.3"
    Pick one species with multiple indexed assemblies (for example two
    *Spiribacter* rows). How does summed abundance across those labels relate to
    the species fraction in `profile.tsv`?

---

## Single-reference mode - `--no-multi`

By default, krepp reports a read against multiple references when it is nearly
equally distant from several. `--no-multi` enforces a single hit:

```bash
krepp dist -i data-new/ref_index -q data-new/query_mixture.fq.gz \
    --no-multi --num-threads 4
```

---

## Distance interpretation

| d_JC | Biological meaning |
|------|---------------------|
| < 0.01 | Same strain or very close ecotype |
| 0.01 - 0.05 | Within-species / ecotype divergence |
| 0.05 - 0.15 | Different species in the same genus |
| 0.15 - 0.30 | Different genus within a family |
| > 0.30 | Family or higher divergence |

Reads from a query whose best reference shares only a **coarse** rank (see
`role` in `query_info.tsv`) often land on that surrogate with distances in the
**within-genus** band or wider, depending on how far apart the taxa are.

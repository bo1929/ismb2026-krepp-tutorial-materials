# Distance Estimation - `krepp dist`

- `krepp dist` computes the Hamming distance between individual sequences (e.g., short reads in our case) in a given FASTA/Q file and sufficiently similar (i.e., <25% similarity) reference genomes in the index.
- These approximate distances are akin to alignment identity: If you were to align this read (e.g., using Bowtie2), the proportion of mismatches are expected to match krepp's distance estimates.
- Individual estimates for reads may have some error, especially at high distances (>15%), but they are accurate on average, and can recover genome-wide similarities.
- These distances are main information krepp relies on for various tasks such as phylogenetic placement, taxonomic classification, and for providing sample-wise (i.e., for the entire file) summaries.

---

## Per-read distances (default)

By default, krepp reports distances between individual sequences and all relevant references:
```bash
krepp dist -i data/toy-index -q data/query_mixture.fq.gz --num-threads 4 -o results/distances_default.tsv
```

The output is in a tab-separated format, consisting of three columns: `SEQ_ID` (sequence/read identifier in the input), `REFERENCE_NAME` (matching reference in the index), `DIST` (estimated distance).

```bash
head results/distances_default.tsv
```
??? question "Expected output:"
    ```
    # software: krepp   version: v0.8.2 invocation :krepp dist -i data/toy-index -q data/query_mixture.fq.gz --num-threads 4 -o results/distances_default.tsv
    SEQ_ID  REFERENCE_NAME  DIST
    NC_013960.1-106607  Nitrosococcus_oceani    0.08852
    NC_013960.1-106607  Nitrosococcus_wardiae   0.08739
    NC_013960.1-106607  Nitrosococcus_watsonii  0.11086
    NC_008818.1-97495   NA NaN
    NC_009663.1-23622   NA NaN
    NC_015559.1-12009   Marinomonas_posidonica  0.08579
    NC_019566.1-5659    NA NaN
    NC_015559.1-3949    Marinomonas_posidonica  0.13930
    ```
**Some observations:**

- The exact order of the rows may change across reads, but the mappings should be identical.
- Notice that, a sequence identifier may appear at multiple consecutive rows, this means that the corresponding read is mapped to multiple references, each with its own distance.
- The number of unique IDs in the first column is equal to the number of reads in the given query file.
- Some rows could have NA; indicating that there is no reference that could be matched for this particular read.

---

## Other options for distance estimation

### Filtering reference hits

By default, krepp reports all reference hits, regardless of their distances.

If you are only interested in references that sufficiently close compared to the best match, use `--filter` flag:
```bash
krepp dist -i data/toy-index -q data/query_mixture.fq.gz --num-threads 4 --filter -o results/distances_filtered.tsv
```
For each read, this only retains mappings that are not too far away from the reference hit with the minimum distance.

??? question "Example impact:"
    The total number of mappings goes down from 80,983 to 59,074 (comparing `grep -v NaN results/distances_*.tsv | wc -l`).
    ```
    # software: krepp   version: v0.8.2 invocation :krepp dist -i data/toy-index -q data/query_mixture.fq.gz --num-threads 4 -o results/distances_default.tsv
    SEQ_ID  REFERENCE_NAME  DIST
    NC_013960.1-106607  Nitrosococcus_oceani    0.08852
    NC_013960.1-106607  Nitrosococcus_wardiae   0.08739
    NC_013960.1-106607  Nitrosococcus_watsonii  0.11086
    NC_008818.1-97495   NA NaN
    NC_009663.1-23622   NA NaN
    NC_015559.1-12009   Marinomonas_posidonica  0.08579
    NC_019566.1-5659    NA NaN
    NC_015559.1-3949    Marinomonas_posidonica  0.13930
    ```

??? question "Q 4.1"
    For a read simulated from a query whose `novelty_level` in `query_info.tsv` is coarse
    (e.g. `kingdom` or `family`), which reference label tends to win `--filter`?
    Does the mash / k-mer distance match your expectation from the phylogeny?

---

### Community abundance - `--summarize`

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
    assembly **B** (another strain, species in the same genus, or a more
    distant best match with a coarser shared rank, e.g. `genus` or `kingdom`), **`krepp dist --summarize`**
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

### Single-reference mode - `--no-multi`

By default, krepp reports a read against multiple references when it is nearly
equally distant from several. `--no-multi` enforces a single hit:

```bash
krepp dist -i data-new/ref_index -q data-new/query_mixture.fq.gz \
    --no-multi --num-threads 4
```

---

## Interpretation of results

| d_JC | Biological meaning |
|------|---------------------|
| < 0.01 | Same strain or very close ecotype |
| 0.01 - 0.05 | Within-species / ecotype divergence |
| 0.05 - 0.15 | Different species in the same genus |
| 0.15 - 0.30 | Different genus within a family |
| > 0.30 | Family or higher divergence |

Reads from a query whose best reference shares only a **coarse** rank (see
`novelty_level` in `query_info.tsv`) often land on that surrogate with distances in the
**within-genus** band or wider, depending on how far apart the taxa are.
# software: krepp	version: v0.8.2	invocation :krepp dist -i data/toy-index -q data/query_mixture.fq.gz --num-threads 2 -o results/distances_default.tsv --dist-max 0.33
SEQ_ID	REFERENCE_NAME	DIST
NZ_CP186025.1-22373	Vibrio_anguillarum	0.00875
NC_007484.1-37612	Nitrosococcus_watsonii	0.07074
NC_007484.1-37612	Nitrosococcus_wardiae	0.20504
NC_007484.1-37612	Nitrosococcus_oceani	0.00009
NC_013960.1-34342	Nitrosococcus_wardiae	0.11605
NC_013960.1-34342	Nitrosococcus_oceani	0.18942
NC_009663.1-23622	NaN	NaN
NC_015559.1-12009	Marinomonas_posidonica	0.08581

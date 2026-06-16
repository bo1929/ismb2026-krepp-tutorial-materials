# Indexing reference genomes

Whether you have a collection of reference genomes or prefer a public database, krepp expects them to be indexed using locality-sensitive hashing.
An index is simply a directory storing certain data structures in binary format, together with a human-readable metadata file.

The krepp command-line interface offers a subcommand called `krepp index` that builds an index from a collection of reference genomes.

!!! danger "Ready-to-use indexes"
    We also offer a [catalogue](https://github.com/bo1929/krepp/wiki/Available-reference-indexes) of ready-to-use indexes constructed from different public databases, available for download.
    In either case, you build or download it once, then perform queries and reuse it for your analyses.

!!! type "Goals in this part of the tutorial"
    - **Build a small index** from a selection of 31 microbial genomes detailed in the previous section
    - **Download an index** built from the Web of Life (WoL) database (>10,000 genomes with a reference phylogeny)

---

## Indexing genomes from scratch

### Input files

Given a set of reference genomes, krepp builds an index via `krepp index -o /path/to/index -i /path/to/input.tsv`.

There are only two required options:

* `-o` (`--index-dir`), which is the output directory path in which the index will be stored,

* `-i` (`--input-file`) a two-column TSV file mapping reference IDs to paths/URLs for FASTA/Q files (gzip compatible).

Each reference (whether a complete assembly or a genome skim) is a file identified by a unique ID.
The input map links these IDs to their file paths using a two-column, tab-separated format: `REF-ID<TAB>/path/to/ref`.
The paths can be either absolute or relative to the current working directory.

For the toy dataset, the input map can be seen below:
```bash
head -5 data/input_map.tsv
```

```
Pyrodictium_occultum	data/reference_genomes/GCF_000007085.1.fna.gz
Thermus_thermophilus	data/reference_genomes/GCF_000006785.2.fna.gz
Allochromatium_vinosum	data/reference_genomes/GCF_900119735.1.fna.gz
Candidatus_Endolissoclinum_sp	data/reference_genomes/GCA_965228445.1.fna.gz
Psychrobacter_cryohalolentis	data/reference_genomes/GCF_000013905.1.fna.gz
```

!!! danger ""
    The IDs are particularly important for phylogenetic placement: the leaves of the phylogeny provided for placement (i.e., a Newick file) must match the IDs given in this file.
    For distance estimation, there is no such requirement and arbitrary IDs can be used.
    <!-- Note that the indexes we make available have reference IDs compatible with the corresponding phylogeny we provide. -->

### Building the index

Run from inside `data` so paths in `input_map.tsv` resolve.
```bash
cd data
krepp index \
    -i input_map.tsv \
    -o toy-index \
    --num-threads 4
cd ..
```
??? question "Expected output:"
    ```
    krepp version: v0.8.2
    Invocation: krepp index -i input_map.tsv -o toy-index --num-threads 4
    Tue May 26 05:20:01 2026
    Reading the tree and initializing the index...
    No tree has given as a guide, the color index could be suboptimal.
    Building the index...
    Internal node: 60	size: 13408544	progress: 61/61
    Finished indexing, elapsed: 19.8957 sec
    Skipped saving a backbone for the index!
    Done converting & saving, elapsed: 0.069542 sec
    Tue May 26 05:20:21 2026
    ```
You will, of course, have different dates and exact elapsed times.

??? more "Multi-threaded indexing"
    The above command took ~27 seconds with a single thread.
    Using multiple threads can give significant speed-ups for large reference sets, but this is optional.
    Note that memory use increases with more threads, so you may want to reduce the thread count if you run out of memory during indexing.
    It does not affect memory use at query time, once indexing is finished.

### Inspecting and verifying an index
We can verify an index by running
```bash
krepp inspect -i data/toy-index | head -n 16
```
??? question "Expected output:"
    ```
    Backbone tree: NA
    ======= Partial index: 0 =======
    krepp version: v0.8.2
    date: 2026-05-26 05:27:50
    seed: 0
    k: 29
    w: 35
    h: 13
    m: 4
    frac: true
    ppos_v: [28, 26, 25, 22, 16, 14, 12, 11, 7, 5, 4, 3, 1]
    npos_v: [0, 2, 6, 8, 9, 10, 13, 15, 17, 18, 19, 20, 21, 23, 24, 27]
    nrows: 33554432
    total_num_kmers: 13408544
    sdust-t: 0
    sdust-w: 0
    ```
This command also reports the parameters used during indexing for future reference.
Alternatively, to display similar information about the index, we can check the metadata file.
```
cat data/toy-index/metadata-*.txt
```

---

### (Optional) Indexing options and parameters

| Flag | Default | Meaning |
|------|---------|---------|
| `-k` (`--kmer-len`) | 29 | Length of *k*-mers [29] |
| `-w` (`--win-len`) | 35 | Length of minimizer window, -w must be greater than or equal to -k [k+6] |
| `-h` (`--num-positions`) | 13 | Number of positions for the LSH [k-16] |

!!! note
    Defaults are calibrated with short reads and microbial genomes in mind, and should perform well for most cases.

If you only want to test different parameters and options without overwriting the existing index, use a different `-o` directory.

??? more "Advanced: batching for ultra-large datasets and FracMinHash"
    Peak memory usage during index construction can also be controlled by partitioning the index into smaller pieces.
    This is done by a FracMinHash-based batching, controlled by options `-m` and `-r`.
    `krepp` partitions the index into `-m` (more or less) equally sized pieces; these partitions can be built independently but queried together.
    The `-r` option determines which partition is constructed: if `--no-frac` is given, only the `-r`th partition is built; otherwise all partitions from 0th to `-r`th are built and saved (see `krepp index --help` for details).
    You don't need to construct all partitions; `krepp` will search in whatever is available, and these partitions can be distributed independently.
    The **default** is `-m 4 -r 1 --frac`, so 50% (0th and 1st out of 4 partitions) of the minimized *k*-mers will be indexed.
    The only requirement is keeping the `-m` value (and of course `-i` and `-t`) fixed across all partitions.
    For instance, one can index 3% of the reference *k*-mers and construct a lightweight index by running:
    ```bash
    krepp index -o /path/to/index -i /path/to/input.tsv -t /path/to/tree.nwk -m 100 -r 2 --no-frac
    krepp index -o /path/to/index -i /path/to/input.tsv -t /path/to/tree.nwk -m 100 -r 1 --no-frac
    krepp index -o /path/to/index -i /path/to/input.tsv -t /path/to/tree.nwk -m 100 -r 0 --no-frac
    ```
    or, alternatively
    ```bash
    krepp index -o /path/to/index -i /path/to/input.tsv -t /path/to/tree.nwk -m 100 -r 2 --frac
    ```
    This is faster but uses more memory during index construction, despite resulting in an index of the same size.
    If this doesn't work well for your task, you can keep adding batches partially and make index larger:
    ```bash
    krepp index -o /path/to/index -i /path/to/input.tsv -t /path/to/tree.nwk -m 100 -r 3 --no-frac
    krepp index -o /path/to/index -i /path/to/input.tsv -t /path/to/tree.nwk -m 100 -r 4 --no-frac
    ```
    or, if you have enough memory, you can do this in parallel:
    ```bash
    seq 3 99 | xargs -I{} -P4 \
        krepp index -o /path/to/index -i /path/to/input.tsv -t /path/to/tree.nwk -m 100 -r {} --no-frac
    ```

??? more "Advanced: using a guide tree during indexing for better compression"
    If you have a phylogeny available (in Newick format), run `krepp index` with the `-t` option and the path to the Newick tree (inside `data/`):
    ```bash
    krepp index \
        -i input_map.tsv \
        -o toy-index \
        -t reference_tree.nwk \
        --num-threads 4
    ```
    An index constructed with a backbone can be used for phylogenetic placement without providing a tree later, and distance estimates are unaffected.
    You can always specify a backbone with the `krepp place` subcommand, using the `-t` (`--nwk-file`) option.
    Note that the tip labels of the phylogeny must match the reference IDs given in `-i /path/to/input.tsv` for `krepp index`.

---

## Installing pre-built indexes

We will use the Web of Life index (tiny version <2 GB) which includes a microbial phylogeny (>10,000 leaves).
!!! tip ""
    You do not need to run the command below if you already followed the steps in the **Setup** section.
    Note that this index is intentionally made extremely lightweight for this tutorial and may underperform for complex samples compared to larger indexes.
```bash
wget --no-check-certificate https://ter-trees.ucsd.edu/data/krepp/index-WoLv1-tiny.tar
mkdir -p data
tar -xf index-WoLv1-tiny.tar -C data/
```
The downloaded directory contains the index.
This particular index is small enough to fit in your laptop's memory.

Run `cat data/index-WoLv1-tiny/metadata-*.txt` to see the configuration and verify the download.
??? question "Expected output:"
    ```
    krepp version: v0.8.2
    date: 2026-05-29 04:53:32
    seed: 0
    k: 27
    w: 33
    h: 12
    m: 36
    frac: false
    ppos_v: [26, 25, 24, 22, 21, 17, 14, 8, 7, 5, 3, 2]
    npos_v: [0, 1, 4, 6, 9, 10, 11, 12, 13, 15, 16, 18, 19, 20, 23]
    nrows: 466034
    total_num_kmers: 184716530
    sdust-t: 0
    sdust-w: 0
    ```

See [this page](https://github.com/bo1929/krepp/wiki/Available-reference-indexes) for the full list of available databases.

??? more "How to choose an index?"
    The available indexes are mostly microbial and the majority come with a reference phylogeny (for placement).
    Some major instances are listed below; they vary in size to fit the memory available to you:

    * **Web of Life (WoL):** lightweight and a reliable phylogeny (67 GB, 41 GB, **2 GB** ← we downloaded this one)

    * **RefSeq-microbial:** denser sampling and an ultra-large [uDance phylogeny](https://doi.org/10.1038/s41587-023-01868-8) (180 GB)

    * **GTDB:** the latest release, an up-to-date reference (230 GB, 100 GB)
    <p>
    <img src="figures/ref-logos.png" alt="References." style="display: block; margin: 1rem auto; max-width: 380px;" />
    </p>
    For any reference database, the rule of thumb is to use the largest index that fits your available memory.
    The indexes vary in size but are based on the same genome set. The level of compression and the extent of *k*-mer downsampling differ.
    Depending on the query data and the task, a small index might perform as accurately as a large one.
    The downloads can be either through our servers or an AWS S3 bucket.


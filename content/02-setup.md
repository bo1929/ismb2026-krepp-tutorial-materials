# Setup

## Hardware requirements

Participants should have access to a machine with at least 8 GB of RAM and 8 GB of free disk space.
You can either run the tutorial on your personal laptop, or use a remote server that you have access to.
Some commands we will run may benefit from a multi-threaded setup, but this is only optional.

---

## Software requirements

We will be using **krepp** version v0.8.2.

| Tool | Version | Purpose |
|------|---------|---------|
| conda / micromamba | any recent | environment management |
| **krepp** | 0.8.2 | distance estimation and placement |
| gappa | 0.9.0 | phylogenetic placement visualization and processing |
| wget / curl | any recent | downloading datasets |

We recommend using **micromamba**, refer to [this link](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) for installation.

---

## Create an environment and install krepp

Set up an environment and install `krepp` (version 0.8.2) via Bioconda, then verify the installation:
```bash
micromamba create -n krepp-tutorial -c conda-forge -c bioconda krepp=0.8.2 -y
micromamba activate krepp-tutorial
krepp --help
```
??? question "Expected output:"
	```
	krepp version: v0.8.2
	krepp: a tool for k-mer-based search, distance estimation & phylogenetic placement.
	Usage: krepp [OPTIONS] SUBCOMMAND

	Options:
	--help
	--verbose,--no-verbose{false}
								Increased verbosity and progress report.
	--seed UINT                 Random seed for the LSH and other parts that require randomness. [0]
	--num-threads UINT          Number of threads to use in OpenMP-based parallelism. [1]

	Subcommands:
	index                       Build an index from k-mers of reference genomes.
	place                       Place queries on a tree with respect to an index.
	dist                        Estimate distances of queries to genomes in an index.
	inspect                     Display statistics and information for a given index.
	sketch                      Create a sketch from k-mers in a single FASTA/FASTQ file.
	seek                        Seek query sequences in a sketch and estimate distances.
	```

Install helper tools into the same environment, and verify:
```bash
micromamba install -c conda-forge -c bioconda gappa=0.9.0 wget -y
gappa --help
```
??? question "Expected output:"
	```
												....      ....
												'' '||.   .||'
													||  ||
													'|.|'
		...'   ....   ... ...  ... ...   ....        .|'|.
		|  ||  '' .||   ||'  ||  ||'  || '' .||      .|'  ||
		|''   .|' ||   ||    |  ||    | .|' ||     .|'|.  ||
		'....  '|..'|'. ||...'   ||...'  '|..'|.    '||'    ||:.
		'....'          ||       ||
					''''     ''''   v0.9.0 (c) 2017-2025
									by Lucas Czech and Pierre Barbera

	Usage: gappa [OPTIONS] SUBCOMMAND

	Options:
	--help FLAG                 Print this help message and exit.
	--version FLAG              Print the gappa version and exit.

	Subcommands:
	analyze                     Commands for analyzing and comparing placement data, that is, finding differences and patterns.
	edit                        Commands for editing and manipulating files like jplace, fasta or newick.
	examine                     Commands for examining, visualizing, and tabulating information in placement data.
	prepare                     Commands for preparing and preprocessing of phylogenetic and placement data.
	simulate                    Commands for random generation of phylogenetic and placement data.
	tools                       Auxiliary commands of gappa.

	gappa - a toolkit for analyzing and visualizing phylogenetic (placement) data
	```

---

## Download and prepare the tutorial data

```bash
git clone https://github.com/bo1929/ismb2026-krepp-tutorial-materials.git
cd ismb2026-krepp-tutorial-materials
bash scripts/setup.sh
```

!!! note
    `setup.sh` will download a small (~2 GB) microbial index.

---

## Verify the inputs are ready

```bash
ls data/
```

??? question "Expected output:"
    ```
    input_map.tsv
    profile.tsv
    profile_species_accessions.tsv
    query_genomes
    query_info.tsv
    query_mixture.fq.gz
    reference_genomes
    reference_info.tsv
    ```

| Path | Description |
|------|-------------|
| `reference_genomes/` | 31 reference genome assemblies (FASTA, gzip-compressed) |
| `reference_info.tsv` | Metadata describing reference labels, taxonomic groups and lineages, and accessions |
| `input_map.tsv` | Two-column TSV mapping reference labels to genome file paths (31 lines) |
| `query_genomes/` | Genome assemblies for the 20 query organisms (FASTA, gzip-compressed) |
| `query_info.tsv` | Metadata for query organisms: taxon and accession |
| `query_mixture.fq.gz` | 100,000 Illumina reads simulated from the query genomes (FASTQ, gzip-compressed) |
| `profile.tsv` | Taxonomic profile of the mock community across taxa (77 lines) |

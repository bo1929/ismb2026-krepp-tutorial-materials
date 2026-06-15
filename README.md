# Data for ISMB 2026 krepp tutorial

Datasets for the hands-on exercises for **metagenomic distance estimation and phylogenetic placement** with [krepp](https://github.com/bo1929/krepp).
This document describes how the repository is organized and what data are provided.
It does not discuss the content of the tutorial.

> [!NOTE]
> This README file is not necessarily helpful for the tutorial, feel free to not read.
> Following the provided instructions and the interactive tutorial will be sufficient.

---

## Repository layout

| Path | Role |
|------|------|
| `data/` | **Tutorial dataset**: profiles, metadata tables, reference/query FASTAs, simulated reads (see [Data overview](#data-overview)). |
| `results/` | Scratch directory for command outputs during the workshop (created by `scripts/setup.sh`). |
| `figures/` | Static images referenced from Markdown (dataset overview plots, indexing logos, etc.). |
| `scripts/` | Setup, data construction, HTML build, and figure generation. |
| `content/` | Tutorial source: one Markdown file per lesson (see [Content structure](#content-structure)). |
| `config.yml` | Site metadata and navigation order for the HTML build. |
| `pages/` | Generated HTML (not committed): `full.html`, per-lesson pages, `index.html`. |
| `SYNTAX.md` | Maintainer notes on Markdown dialect and admonitions. |
| `BUILD.md` | How to build the HTML site. |

Large or generated assets are listed in `.gitignore` (indexes built during the course, `pages/`, local `results/`, etc.).

---

## Content structure

All lessons live under `content/` as numbered Markdown files. Order and sidebar grouping are defined in `config.yml` under `nav` (currently a single block, **Tutorial**).

| File | Lesson (nav title) | Role in the flow |
|------|-------------------|------------------|
| `01-overview.md` | Overview | Concepts and overview of the method. |
| `02-setup.md` | Setup | Environment, krepp install, downloading the WoL tiny index, expected `data/` layout. |
| `03-dataset.md` | Dataset | Describes the mock community and points to every `data/` artifact. |
| `04-indexing.md` | Indexing | WoL index usage and building the small toy index from `input_map.tsv`. |
| `05-distance-estimation.md` | Distance Estimation | `krepp dist` on the query mixture against an index. |
| `06-phylogenetic-placement.md` | Phylogenetic Placement | `krepp place` on the same reads and index. |
| `07-taxonomic-placement.md` | Taxonomic Placement | Taxonomic mapping and placement. |
| `08-visualization.md` | Visualization | Inspecting and visualizing results. |
| `09-conclusions.md` | Conclusions | Wrap-up and summary. |

Authoring conventions (admonitions, code fences, tables) are documented in `SYNTAX.md`.

### Building the HTML site

See `BUILD.md`.

## Data overview
**Mock metagenome with known truth:** 20 profile species (marine Bacteria/Archaea), 100k simulated reads (`query_mixture.fq.gz`), 20 query assemblies, 31 reference genomes.
Some reads come from query genomes that are novel, in other words, without a close representation among the small reference genome set.

| What | File(s) in `data/` |
|------|-----------|
| Ground truth abundances | `profile.tsv` |
| Metadata for the query genomes | `query_info.tsv`, `query_taxonomy.tsv` |
| Query genomes \& simulated reads | `query_genomes/`, `query_mixture.fq.gz` |
| Metadata for the reference genomes | `reference_info.tsv`, `reference_taxonomy.tsv` |
| Reference genomes, their IDs, and the reference phylogeny | `input_map.tsv`, `reference_genomes/`, and `reference_tree.nwk` |
| WoL-v1 tiny index | `bash scripts/setup.sh` downloads `data/index-WoLv1-tiny/` (~10,000 refs.) |

**Quick start:** `bash scripts/setup.sh` then `bash scripts/build_site.sh` and open `pages/index.html`.

Install krepp, following `content/02-setup.md`.

Build and/or download indexes following `content/04-indexing.md` (also see [krepp wiki](https://github.com/bo1929/krepp/wiki/Available-reference-indexes)).

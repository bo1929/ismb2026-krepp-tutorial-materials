# krepp Tutorial: Marine Metagenome Placement

**ISMB 2026** | **krepp v0.8.0** | **Ali Osman Berk Sapci**

Simulated marine sample: mixed references plus *Prochlorococcus* AS9601 held out to learn **dist**, **placement**, **taxonomy**, and **gappa** heat maps.

```
genomes --> [krepp index] --> index ---> [dist --summarize] --> abundance
                                       \-> [place] -----------> placements / taxonomy
```

| # | Topic | File |
|---|-------|------|
| 1-10 | Linear workshop | `content/01-overview.md` through `10-cheatsheet.md` |
| + | Cases B/C/D | [`content/cases/`](content/cases/) |
| key | Instructor answers | [`content/answers.md`](content/answers.md) |

## Quick commands

```bash
git clone <repo>; cd ismb2026-krepp-tutorial
bash scripts/setup.sh
cd data && krepp index -i input_map.tsv -o index -t tree.nwk --num-threads 4 && cd ..
mkdir -p results
krepp dist  -i data/index -q data/query.fq --summarize --num-threads 4
krepp place -i data/index -q data/query.fq --summarize --num-threads 4
```

Large assets (`data/genomes/`, `index/`, reads) stay **gitignored** unless you extract Zenodo payloads (URLs TBD) or rerun `scripts/setup.sh`. Use `precomputed/` when you skip krepp entirely.

## Repo layout

- `content/` markdown (see [`SYNTAX.md`](SYNTAX.md) for callout syntax).
- `config.yml` nav for `scripts/build_html.py`; `tutorial.html` = `bash scripts/build_site.sh`.
- `scripts/` orchestration (`setup.sh`, genome fetch / tree / simulator / figures / HTML generator).
- `precomputed/` saved outputs cases base + B/C/D.

## Cases (beyond base worksheet)

| Label | Twist |
|-------|-------|
| B | Novel at 0.5% abundance |
| C | Two simultaneous holdouts |
| D | Single-flank strain holdout |

## Build handbook

```bash
pip install pymdown-extensions pyyaml markdown
bash scripts/build_site.sh
```

Outputs `tutorial.html` (gitignored) from `content/` + `scripts/template.html`.

## Reference

Sapci & Mirarab, *Genome Biology* 27:108 (2026).
https://doi.org/10.1186/s13059-026-03999-y

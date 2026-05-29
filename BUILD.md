## 
## Repo layout

- `content/` markdown (see [`SYNTAX.md`](SYNTAX.md) for callout syntax).
- `config.yml` nav for `scripts/build_html.py`; `tutorial.html` = `bash scripts/build_site.sh`.
- `scripts/` orchestration (`setup.sh`, genome fetch / tree / simulator / figures / HTML generator).
- `precomputed/` saved outputs cases base + B/C/D.

## Build handbook

```bash
pip install pymdown-extensions pyyaml markdown
bash scripts/build_site.sh
```

Outputs `tutorial.html` (gitignored) from `content/` + `scripts/template.html`.

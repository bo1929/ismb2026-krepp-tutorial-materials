# Repo layout

- `content/` markdown (see [`SYNTAX.md`](SYNTAX.md) for callout syntax).
- `config.yml` nav for `scripts/build_html.py`; `bash scripts/build_site.sh` writes `pages/`.

# Build handbook

```bash
pip install pymdown-extensions pyyaml markdown
bash scripts/build_site.sh
```

Outputs under `pages/` (gitignored):

- `full.html`: full standalone page with sidebar (`scripts/template.html`).
- `<stem>.html`: isolated pages per section, no sidebar (`scripts/page_template.html`).
- `index.html`: index for sections (`scripts/index_template.html`).

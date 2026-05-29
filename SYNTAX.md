# Tutorial Markdown syntax and HTML mapping

[!CAUTION]
> **Participants should ignore this file.**
> This document is **for maintainers and authors** who edit files under `content/`.

---

## Pipeline

Tutorial pages are Markdown files under `content/`.
The handbook and all the pages are built by `scripts/build_html.py`, which:

1. Reads `config.yml`: `docs_dir`, `nav`.
Preferred shape is **grouped lists** (each top item is `SectionName:` with indented `- Title: relative/path.md`).
`nav_badges` scopes badges by that section heading, then YAML title keys (reuse names like `Overview`).
Older styles (`level:` rows, flat `section:`, or fully flat lists) still parse when no grouped lists appear.
Optional `handbook`, legacy `sidebar_nav`, and `site_*`.

2. Loads the HTML shell from `scripts/template.html`.
3. Converts each Markdown file with Python-Markdown plus PyMdown extensions.

The dialect is whatever those extensions define, not arbitrary CommonMark/GitHub flavored Markdown alone.

---

## Extensions loaded (`scripts/build_html.py`)

| Extension | Role |
|-----------|------|
| `tables` | GitHub-style pipe tables. |
| `fenced_code` | Fenced triple-backtick code blocks (also used with superfences). |
| `admonition` | `!!!` callout blocks. |
| `markdown.extensions.toc.TocExtension(permalink=False)` | Table of contents / heading anchors; permalink symbols next to headings are **off**. |
| `pymdownx.details` | `???` collapsible sections. |
| `pymdownx.superfences` | Enhanced fenced blocks (pairs with highlighting). |
| `pymdownx.highlight` with `use_pygments: False` | Non-Pygments code output used by our post-processing. |

**Not enabled:** attribute lists (`{: ...}`) and other Markdown-Extra shortcuts are unavailable unless you paste raw HTML (not documented as supported for authors).

**Upstream docs:** Python-Markdown admonition: https://python-markdown.github.io/extensions/admonition/

PyMdown Details extension: https://facelessuser.github.io/pymdown-extensions/extensions/details/

---

## Standard Markdown used in `content/`

- **Headings:** `#` ... `######` (ATX style).
- **Emphasis:** `**bold**`, `*italic*`.
- **Inline code:** `` `like this` ``.
- **Links:** `[label](relative-or-absolute-url)`.
- **Lists:** unordered (`-`, `*`) and ordered (`1.`).
- **Horizontal rule:** `---` on its own line.
- **Paragraphs:** blank line between blocks.

---

## Tables

Authors use GitHub-flavored pipe tables:

```markdown
| Col A | Col B |
|-------|-------|
| one   | two   |
```

**HTML mapping:** Each `<table>...</table>` is wrapped in `<div class="table-wrap">...</div>` by `scripts/build_html.py` (`wrap_tables()`) so wide tables scroll in the handbook layout.

---

## Fenced code blocks

````markdown
```bash
echo hello
```
````

Without a language tag, the fence still produces a `<pre>/<code>` block.

With a language identifier (here `bash`), PyMdown produces `<pre class="highlight"><code class="language-bash">...</code></pre>`.

**Post-processing (`add_data_lang()`):** opening tags are rewritten to `<pre data-lang="bash"><code>...</code></pre>` (or `<pre><code>...</code></pre>` with no lang).
**`scripts/template.html`** uses `pre[data-lang="..."]` for the visible language badge in the handbook.

---

## Admonitions: `!!!`

Syntax:

```markdown
!!! type "Optional title shown in uppercase"
    Every line of the body is indented (four spaces is the usual convention).
    Blank lines inside the body are allowed if indented.
```

Lines after `!!!` that belong to the admonition must stay **indented** until the block ends.

**Title omitted:** valid for some types:

```markdown
!!! note
    Body only. Python-Markdown still emits a title line; for `note` it is typically the word **Note**.
```

Example in this repo: `content/04-indexing.md`.

### Types used under `content/`

These are the `type` values:

`abstract`, `info`, `note`, `tip`, `success`, `warning`, `danger`, `next`.

For instance, use `!!! next "Your title"` (or bare `!!! next` with a default **Next** title) to mark the next step or highlight the immediate next tutorial action so it scans like an anchor among normal prose.

### HTML emitted

Rough structure:

```html
<div class="admonition {type}">
  <p class="admonition-title">Title Here</p>
  <p>...</p>
</div>
```

The second element is omitted or adjusted if there is only body text; trust the converter output when debugging.

### Styling (`scripts/template.html`)

Base styling: `.admonition`, `.admonition-title`.

Variant grouping (same border and background cues):

| CSS selector | Typical use in content |
|----------------|-------------------------|
| `.admonition.note`, `.admonition.info` | Neutral / informational |
| `.admonition.tip`, `.admonition.success` | Positive / teaching tips |
| `.admonition.warning` | Caution |
| `.admonition.danger`, `.admonition.abstract` | Strong callouts and abstract-style intros |
| `.admonition.next` | Next step / resume-here cue (accent stripe and ring-shaped marker in `template.html`) |

Untyped visuals fall through to the generic `.admonition` blues.

---

## Collapsible blocks: `???` (PyMdown Details)

### Questions (exercises)

Tutorial exercises use the `question` class:

```markdown
??? question "Q 2.1"
    Prompt or answer text, indented under the headline.
```

- **Collapsed by default.** Readers click the summary to expand.
- The keyword `question` becomes part of emitted classes (`question` alongside optional `admonition` markup depending on PyMdown version/settings).

Rough HTML:

```html
<details class="question">
  <summary>Q 2.1</summary>
  <p>...</p>
</details>
```

`scripts/template.html` styles `details.question` and `details.admonition.question` with a green `?` badge on the summary line.

### Optional detail (`more`)

Use **three question marks** `???`, not `!!!`.
The latter is an **admonition** and is **not** collapsible; optional blocks must use **Details**:

For **extra background** or details that should stay out of the main path (like a spoiler or deep dive), use `more`:

```markdown
??? more "Why this matters"
    Longer explanation, optional links, assumptions, etc.
```

Rough HTML:

```html
<details class="more">
  <summary>Why this matters</summary>
  <p>...</p>
</details>
```

`scripts/template.html` styles `details.more` with an accent-tone `+` / `-` control (collapsed vs open), distinct from quiz-style question blocks.

---

## Author checklist

1. Indent admonition and `???` (`question`, `more`, etc.) bodies consistently (four spaces is standard).
2. Leave a blank line before and after fenced code and callouts where Markdown looks ambiguous.
3. Sync `config.yml` `nav`: group pages under `SectionName:` blocks; keep `nav_badges` aligned with those headings and link titles.
4. Rebuild locally: `bash scripts/build_site.sh` or `python3 scripts/build_html.py --root .`

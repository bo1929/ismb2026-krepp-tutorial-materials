#!/usr/bin/env python3
"""
build_html.py - merge content/*.md (per config.yml nav) into tutorial.html.
Uses python-markdown + pymdown; post-processes fenced <pre>: always emits data-lang
(empty string when no fence language); table-wrap.
"""

import re
import sys
import json
import html as html_escape
import argparse
from pathlib import Path

import yaml
import markdown
from markdown.extensions.toc import TocExtension

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

EXTENSIONS = [
    "tables",
    "fenced_code",
    "admonition",
    TocExtension(permalink=False),
    "pymdownx.details",
    "pymdownx.superfences",
    "pymdownx.highlight",
]

EXTENSION_CONFIGS = {
    # Use_pygments=False gives <code class="language-X"> which we can pick up
    # with data-lang post-processing; avoids a pygments dependency at runtime.
    "pymdownx.highlight": {
        "use_pygments": False,
    },
    "pymdownx.superfences": {
        "disable_indented_code_blocks": False,
    },
}


def _nav_uses_level_schema(nav):
    """True when every ``nav`` item is a dict with a numeric ``level`` (section rows use ``section:`` inside)."""
    if not isinstance(nav, list) or not nav:
        return False
    for entry in nav:
        if not isinstance(entry, dict) or "level" not in entry:
            return False
    return True


def _nav_has_group_lists(nav):
    """True when any item is ``{ "Section": [ nested ... ] }``."""

    if not isinstance(nav, list):
        return False
    for e in nav:
        if isinstance(e, dict) and len(e) == 1:
            v = next(iter(e.values()))
            if isinstance(v, list):
                return True
    return False


def flatten_nav_level(nav, docs_dir):
    """Level-based YAML: ``- level: N`` plus ``section: heading`` OR one ``NavTitle: relative.md``.

    Subsequent pages reuse the nearest preceding ``section`` name for scoped ``nav_badges``.
    """
    items = []
    current_section = ""

    if not isinstance(nav, list):
        return items

    for entry in nav:
        if not isinstance(entry, dict):
            continue
        try:
            raw_level = entry["level"]
            nav_level = int(raw_level)
        except (KeyError, TypeError, ValueError):
            continue

        if "section" in entry:
            lab = entry["section"]
            s = "" if lab is None else str(lab).strip()
            if s:
                items.append((s, None, "", nav_level))
                current_section = s
            continue

        page_keys = [k for k in entry.keys() if k not in ("level", "section")]
        if len(page_keys) != 1:
            continue
        pk = page_keys[0]
        pv = entry[pk]
        if not isinstance(pv, str):
            continue
        rel = pv.strip()
        items.append((pk, docs_dir / rel, current_section, nav_level))

    return items


# Legacy: single-key ``{"section": "..."}`` row (divider) vs ``Title: *.md``.
NAV_SECTION_KEY = "section"


def _nav_value_looks_like_md_relpath(val):
    v = val.strip()
    if not v:
        return False
    low = v.lower()
    if low.endswith((".md", ".markdown")):
        return True
    if "/" in v or "\\" in v:
        return True
    return False


def _walk_hierarchical_children(children, docs_dir, items, badge_scope, group_depth):
    """Recurse subsection lists; ``badge_scope`` is the YAML section used for scoped ``nav_badges``."""

    if not isinstance(children, list):
        return
    page_level = group_depth + 1
    for child in children:
        if not isinstance(child, dict) or len(child) != 1:
            continue
        ck, cv = next(iter(child.items()))
        if isinstance(cv, str):
            rel = cv.strip()
            if rel:
                items.append((ck, docs_dir / rel, badge_scope, page_level))
        elif isinstance(cv, list):
            sub = ck.strip()
            if sub:
                gdepth = group_depth + 1
                items.append((sub, None, "", gdepth))
                _walk_hierarchical_children(cv, docs_dir, items, sub, gdepth)


def flatten_nav_hierarchical(nav, docs_dir):
    """Grouped nav: ``- Core Tutorial:\\n      - Overview: doc.md``

    Lone top-level ``Title: path`` rows attach to the most recent group heading.

    Still accepts legacy plain ``section: Heading`` divider rows alongside groups.
    """
    items = []
    last_section = [""]

    if not isinstance(nav, list):
        return items

    for entry in nav:
        if isinstance(entry, str):
            items.append((entry.strip(), None, "", 0))
            continue

        if not isinstance(entry, dict) or len(entry) != 1:
            continue

        key, value = next(iter(entry.items()))

        if isinstance(value, list):
            heading = key.strip()
            if heading:
                items.append((heading, None, "", 1))
                last_section[0] = heading
                _walk_hierarchical_children(value, docs_dir, items, heading, 1)
            continue

        if isinstance(value, str):
            vs = value.strip()
            if key == NAV_SECTION_KEY:
                if vs and _nav_value_looks_like_md_relpath(vs):
                    items.append((key, docs_dir / vs, last_section[0], 2))
                elif vs:
                    items.append((vs, None, "", 1))
                    last_section[0] = vs
                continue
            if vs:
                items.append((key, docs_dir / vs, last_section[0], 2))

    return items


def _flatten_nav_legacy_list(nav, docs_dir, section_tag):
    """``section_tag`` is a length-1 list: current sidebar section heading for badges."""
    items = []
    if not isinstance(nav, list):
        return items

    for entry in nav:
        if isinstance(entry, str):
            items.append((entry.strip(), None, "", 0))
            continue

        if not isinstance(entry, dict):
            continue

        single = len(entry) == 1
        if single:
            title, value = next(iter(entry.items()))
            if title == NAV_SECTION_KEY and isinstance(value, str):
                vs = value.strip()
                if vs and _nav_value_looks_like_md_relpath(vs):
                    items.append((title, docs_dir / vs, section_tag[0], 0))
                elif vs:
                    items.append((vs, None, "", 0))
                    section_tag[0] = vs
                continue
            if isinstance(value, str):
                items.append((title, docs_dir / value.strip(), section_tag[0], 0))
                continue

            if isinstance(value, list):
                h = title.strip()
                items.append((h, None, "", 0))
                section_tag[0] = h
                items.extend(_flatten_nav_legacy_list(value, docs_dir, section_tag))
                continue

    return items


def flatten_nav_legacy(nav, docs_dir):
    """Legacy YAML: optional ``section:`` divider, ``Title: path``, nested subsection lists."""
    return _flatten_nav_legacy_list(nav, docs_dir, [""])


def flatten_nav(nav, docs_dir):
    if _nav_uses_level_schema(nav):
        return flatten_nav_level(nav, docs_dir)
    if _nav_has_group_lists(nav):
        return flatten_nav_hierarchical(nav, docs_dir)
    return flatten_nav_legacy(nav, docs_dir)


def _nested_nav_badge_roots(nb):
    """True when ``nav_badges`` maps each section heading to an inner badge dict."""
    if not isinstance(nb, dict) or not nb:
        return False
    return all(isinstance(v, dict) for v in nb.values())


def resolve_nav_badge(config, section_heading, nav_title, stem, prefix, idx):
    """Sidebar `.num`: scoped ``nav_badges`` by section then nav YAML key; else stem map; ``sidebar_nav``; defaults."""
    nh = section_heading.strip() if section_heading else ""
    nt = nav_title if nav_title is not None else ""

    def pick(v):
        if v is None:
            return None
        s = str(v).strip()
        return s if s else None

    nb = config.get("nav_badges")
    nested_mode = isinstance(nb, dict) and nb and _nested_nav_badge_roots(nb)

    got = None
    if nested_mode and nh:
        inner = nb.get(nh)
        if isinstance(inner, dict):
            got = pick(inner.get(nt))

    if got is None and isinstance(nb, dict) and not nested_mode:
        for key in (stem, prefix):
            if key is None:
                continue
            got = pick(nb.get(key))
            if got is not None:
                break

    sn = config.get("sidebar_nav")
    if got is None and isinstance(sn, dict):
        lbl = sn.get("badge_labels")
        if isinstance(lbl, dict):
            for key in (stem, prefix):
                if key is None:
                    continue
                gv = pick(lbl.get(key))
                if gv is not None:
                    got = gv
                    break

    if got is not None:
        return got

    if prefix.isdigit():
        return str(int(prefix))
    return str(idx)


def build_nav_html(flat_nav, page_indices, config):
    """Build the inner HTML for <nav id="nav">.

    flat_nav rows: ``(sidebar_title, path_or_None, badge_section_heading, nav_level)``.
    """
    lines = []
    last_heading_emit = None

    for row in flat_nav:
        title, path, badge_section, nav_level = row
        lev_g = nav_level if nav_level > 0 else ""
        lev_a = nav_level if nav_level > 1 else ""

        if path is None:
            display = str(title).strip()
            if display and display != last_heading_emit:
                last_heading_emit = display
                dattr = f' data-level="{lev_g}"' if lev_g else ""
                esc = html_escape.escape(display, quote=False)
                lines.append(f'<div class="nav-group"{dattr}>{esc}</div>')
            continue

        idx = page_indices.get(str(path))
        if idx is None:
            continue

        stem = path.stem
        prefix = stem.split("-")[0] if "-" in stem else stem
        num = resolve_nav_badge(
            config, badge_section or "", title, stem, prefix, idx,
        )

        dattr_a = f' data-nav-level="{lev_a}"' if lev_a else ""
        lines.append(
            f'<a data-page="{idx}"{dattr_a}>'
            f'<span class="num">{num}</span> {title}'
            "</a>"
        )

    return "\n      ".join(lines)



# ---------------------------------------------------------------------------
# Markdown conversion
# ---------------------------------------------------------------------------

def convert_md(md_path):
    """Read a markdown file and return its HTML body."""
    text = Path(md_path).read_text(encoding="utf-8")
    md = markdown.Markdown(
        extensions=EXTENSIONS,
        extension_configs=EXTENSION_CONFIGS,
    )
    return md.convert(text)


# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------

def wrap_tables(html):
    """Wrap bare <table> elements in <div class="table-wrap">."""
    html = re.sub(r"(?<!table-wrap\">)<table", '<div class="table-wrap"><table', html)
    html = re.sub(r"</table>", "</table></div>", html)
    return html


def add_data_lang(html):
    """Convert pymdown fenced <pre class="highlight"> to <pre data-lang="..."><code>.

    Every fenced block carries ``data-lang``, including ``data-lang=""`` when the
    writer omits the fence language, so the handbook shell can rely on uniform
    markup and CSS to hide empty labels."""

    def _replace(m):
        raw = m.group(1)
        body = m.group(2)
        dl = raw.strip() if raw else ""
        if dl:
            aq = html_escape.escape(dl, quote=True)
            return f'<pre data-lang="{aq}"><code>{body}</code></pre>'
        return f'<pre data-lang=""><code>{body}</code></pre>'

    # pymdownx.highlight with use_pygments=False produces:
    #   <pre class="highlight"><code class="language-bash">...</code></pre>
    # or (no language): <pre class="highlight"><code>...</code></pre>
    html = re.sub(
        r'<pre class="highlight"><code(?:\s+class="language-([^"]*)")?>(.*?)</code></pre>',
        _replace,
        html,
        flags=re.DOTALL,
    )
    return html


def ensure_pre_data_lang_attr(html):
    """Add ``data-lang=""`` to bare ``<pre>`` tags missing ``data-lang``."""
    def _inj(m):
        inner = m.group(1).strip()
        if re.search(r"\bdata-lang\s*=", inner, re.I):
            return m.group(0)
        rest = inner
        sep = ""
        if rest:
            sep = " "
        return f'<pre data-lang=""{sep}{rest}>'

    return re.sub(r"<pre([^>]*)>", _inj, html)


def wrap_code_lines(html):
    """Wrap each line inside <pre><code> in <span class="l0/l1"> for zebra striping.

    Per-line spans keep stripe colours correct even when a long line soft-wraps."""

    def _wrap(m):
        pre_attrs = m.group(1)
        code_content = m.group(2)

        if not code_content:
            return m.group(0)

        content = code_content.replace("\r\n", "\n").replace("\r", "\n")
        lines = content.rstrip("\n").split("\n")
        if not lines:
            return m.group(0)

        wrapped_lines = []
        for i, line in enumerate(lines):
            cls = "l0" if i % 2 == 0 else "l1"
            wrapped_lines.append(f'<span class="{cls}">{line}</span>')

        return f"<pre{pre_attrs}><code>{''.join(wrapped_lines)}</code></pre>"

    return re.sub(
        r"<pre([^>]*)><code>(.*?)</code></pre>",
        _wrap,
        html,
        flags=re.DOTALL,
    )


def fix_highlight_pre(html):
    """No-op: highlight div handling is now in add_data_lang."""
    return html


def postprocess(html):
    html = fix_highlight_pre(html)
    html = add_data_lang(html)
    html = ensure_pre_data_lang_attr(html)
    html = wrap_code_lines(html)
    html = wrap_tables(html)
    return html


# ---------------------------------------------------------------------------
# Page HTML builder
# ---------------------------------------------------------------------------

def build_page_html(idx, body_html):
    return (
        f'<div class="page" id="page-{idx}">\n'
        f'<div class="page-inner">\n'
        f'{body_html}\n'
        f'</div>\n'
        f'</div>'
    )


def fill_handbook_placeholders(template, config):
    """Inject handbook shell strings from config. Optional block: missing keys are empty."""
    hb = config.get("handbook")
    if not isinstance(hb, dict):
        hb = {}

    def hb_str(key):
        v = hb.get(key)
        if v is None:
            return ""
        return str(v).strip()

    html_title = hb_str("html_title")
    meta_desc = hb_str("meta_description")
    meta_author = hb_str("meta_author")
    html_lang = hb_str("html_lang")
    sidebar_kicker = hb_str("sidebar_kicker")
    sidebar_heading = hb_str("sidebar_heading")
    sidebar_tagline = hb_str("sidebar_tagline")
    footer_prev = hb_str("footer_prev")
    footer_next = hb_str("footer_next")

    def tx(s):
        return html_escape.escape(s, quote=False)

    repl = {
        "{{HTML_LANG}}": tx(html_lang),
        "{{HTML_TITLE}}": tx(html_title),
        "{{META_DESCRIPTION}}": html_escape.escape(meta_desc, quote=True),
        "{{META_AUTHOR}}": html_escape.escape(meta_author, quote=True),
        "{{SIDEBAR_KICKER}}": tx(sidebar_kicker),
        "{{SIDEBAR_HEADING}}": tx(sidebar_heading),
        "{{SIDEBAR_TAGLINE}}": tx(sidebar_tagline),
        "{{FOOTER_PREV}}": tx(footer_prev),
        "{{FOOTER_NEXT}}": tx(footer_next),
    }

    footer_js = (
        f"const HB_FOOTER_PREV = {json.dumps(footer_prev)};\n"
        f"const HB_FOOTER_NEXT = {json.dumps(footer_next)};"
    )

    for key, val in repl.items():
        template = template.replace(key, val)
    template = template.replace("{{HB_FOOTER_VARS}}", footer_js)
    return template


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=None,
                    help="Root directory (default: parent of scripts/)")
    ap.add_argument("--config", default=None,
                    help="tutorial config YAML (default: <root>/config.yml)")
    args = ap.parse_args()

    scripts_dir = Path(__file__).parent.resolve()
    root = Path(args.root).resolve() if args.root else scripts_dir.parent

    config_path = Path(args.config).resolve() if args.config else root / "config.yml"
    template_path = scripts_dir / "template.html"
    out_path = root / "tutorial.html"

    if not config_path.exists():
        sys.exit(f"config not found: {config_path}")
    if not template_path.exists():
        sys.exit(f"template.html not found at {template_path}")

    cfg_dir = config_path.parent
    with open(config_path) as fh:
        config = yaml.safe_load(fh)

    docs_relative = Path(config.get("docs_dir", "docs"))
    docs_dir = docs_relative if docs_relative.is_absolute() else cfg_dir / docs_relative

    # Flatten nav
    flat_nav = flatten_nav(config.get("nav", []), docs_dir)

    # Assign page indices (only entries with a real path get an index)
    page_entries = [(t, p, sec) for (t, p, sec, _) in flat_nav if p is not None]
    page_indices = {str(p): i for i, (_, p, _) in enumerate(page_entries)}

    # Convert markdown files
    print(f"Building {len(page_entries)} pages from {docs_dir}...")
    pages_html_parts = []

    for idx, (title, md_path, _sec) in enumerate(page_entries):
        if not md_path.exists():
            print(f"  WARNING: {md_path} not found, skipping")
            continue
        body = convert_md(md_path)
        body = postprocess(body)
        pages_html_parts.append(build_page_html(idx, body))
        print(f"  [{idx}] {title}")

    # Build nav
    nav_html = build_nav_html(flat_nav, page_indices, config)

    # Fill template
    template = template_path.read_text(encoding="utf-8")
    template = fill_handbook_placeholders(template, config)
    output = (
        template
        .replace("{{NAV}}", nav_html)
        .replace("{{PAGES}}", "\n\n    ".join(pages_html_parts))
    )

    out_path.write_text(output, encoding="utf-8")
    size_kb = out_path.stat().st_size // 1024
    print(f"\nWrote {out_path}  ({size_kb} KB)")


if __name__ == "__main__":
    main()

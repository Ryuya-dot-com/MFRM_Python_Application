"""Build the static GitHub Pages site from canonical repository documentation.

Uses only the Python standard library so both local previews and GitHub Actions
builds run deterministically with zero extra package installations.
"""
from __future__ import annotations

import argparse
import html
from pathlib import Path
import re
import shutil


GITHUB_BLOB_BASE = "https://github.com/Ryuya-dot-com/MFRM_Python_Application/blob/main/"

PAGE_SPECS = [
    ("docs/desktop_user_guide.en.md", "desktop-guide.html", "en", "MFRM Desktop User Guide — Windows & macOS"),
    ("docs/desktop_user_guide.ja.md", "desktop-guide-ja.html", "ja", "MFRM デスクトップ利用ガイド"),
    ("docs/analysis_scope.md", "analysis-scope.html", "en", "Analyses, Interpretation & Reporting Limits — MFRM"),
    ("docs/apa_manuscript_template.md", "manuscript-guide.html", "en", "APA & Manuscript Reporting Guide — MFRM"),
    ("README.ja.md", "overview-ja.html", "ja", "MFRM — 日本語概要"),
]

ROUTE_MAP = {
    "README.md": "index.html",
    "README.ja.md": "overview-ja.html",
    "docs/desktop_user_guide.en.md": "desktop-guide.html",
    "docs/desktop_user_guide.ja.md": "desktop-guide-ja.html",
    "docs/analysis_scope.md": "analysis-scope.html",
    "docs/apa_manuscript_template.md": "manuscript-guide.html",
    "examples/desktop_writing.csv": "examples/desktop_writing.csv",
}


def github_slug(text: str) -> str:
    clean = re.sub(r"<[^>]+>", "", text)
    clean = re.sub(r"`([^`]*)`", r"\1", clean)
    clean = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", clean)
    clean = re.sub(r"[*_~]", "", clean).strip().lower()
    clean = re.sub(r"[^\w\s-]", "", clean)
    return re.sub(r"\s+", "-", clean)


def resolve_href(raw_url: str, source_md: Path, repo_root: Path) -> str:
    url = raw_url.strip()
    if not url or re.match(r"^[a-zA-Z][a-zA-Z0-9+.-]*:|^//|^mailto:", url):
        return url
    if url.startswith("#"):
        return url
    path_part, sep, frag = url.partition("#")
    target = (source_md.parent / path_part).resolve()
    try:
        rel = target.relative_to(repo_root).as_posix()
    except ValueError:
        return url
    suffix = f"#{frag}" if sep else ""
    if rel in ROUTE_MAP:
        return ROUTE_MAP[rel] + suffix
    if rel.startswith("docs/images/"):
        return "assets/images/" + rel.removeprefix("docs/images/") + suffix
    return GITHUB_BLOB_BASE + rel + suffix


def render_inline(text: str, source_md: Path, repo_root: Path) -> str:
    tokens: list[str] = []

    def stash(fragment: str) -> str:
        idx = len(tokens)
        tokens.append(fragment)
        return f"\x00TOK{idx}\x00"

    # 1. Preserve inline code spans first
    def repl_code(m: re.Match[str]) -> str:
        return stash(f"<code>{html.escape(m.group(1), quote=False)}</code>")

    text = re.sub(r"`([^`]+)`", repl_code, text)

    # 2. Preserve allowed inline HTML tags (<a ...>, <img ...>, <br ...>, <strong>, <em>, <code>)
    def repl_html_tag(m: re.Match[str]) -> str:
        tag = m.group(0)
        tag = re.sub(
            r'(\b(?:href|src)=")([^"]+)(")',
            lambda am: am.group(1) + resolve_href(html.unescape(am.group(2)), source_md, repo_root) + am.group(3),
            tag,
        )
        return stash(tag)

    text = re.sub(r"</?(?:a|img|br|strong|em|b|i|span|sup|sub)\b[^>]*>", repl_html_tag, text)

    # 3. Images ![alt](src)
    def repl_img(m: re.Match[str]) -> str:
        alt = html.escape(m.group(1), quote=True)
        src = html.escape(resolve_href(m.group(2), source_md, repo_root), quote=True)
        return stash(f'<img src="{src}" alt="{alt}" loading="lazy" />')

    text = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)", repl_img, text)

    # 4. Links [label](url)
    def repl_link(m: re.Match[str]) -> str:
        label = render_inline(m.group(1), source_md, repo_root)
        href = html.escape(resolve_href(m.group(2), source_md, repo_root), quote=True)
        ext = ' target="_blank" rel="noopener"' if href.startswith("http") else ""
        return stash(f'<a href="{href}"{ext}>{label}</a>')

    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", repl_link, text)

    # Escape remaining plain text
    escaped = html.escape(text, quote=False)

    # Bold & Italic
    escaped = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", escaped)
    escaped = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<em>\1</em>", escaped)

    for idx, frag in enumerate(tokens):
        escaped = escaped.replace(f"\x00TOK{idx}\x00", frag)
    return escaped


def split_table_row(line: str) -> list[str]:
    s = line.strip()
    if s.startswith("|"):
        s = s[1:]
    if s.endswith("|"):
        s = s[:-1]
    return [cell.strip() for cell in s.split("|")]


def is_table_divider(line: str) -> bool:
    s = line.strip()
    if "|" not in s:
        return False
    cells = split_table_row(s)
    return bool(cells) and all(bool(re.fullmatch(r":?-{3,}:?", c)) for c in cells)


def render_markdown(md_text: str, source_md: Path, repo_root: Path) -> tuple[str, list[tuple[int, str, str]]]:
    lines = md_text.splitlines()
    out: list[str] = []
    toc: list[tuple[int, str, str]] = []
    slug_counts: dict[str, int] = {}
    i = 0
    n = len(lines)

    while i < n:
        line = lines[i]
        stripped = line.strip()

        if not stripped:
            i += 1
            continue

        # Raw block HTML (<p ...>, <h1 ...>, <a id=...>, <details>, </details>, <summary>)
        if re.match(r"^<(?:p|h[1-6]|div|details|summary)\b", stripped) or re.match(r'^<a\s+id="[^"]+"\s*></a>$', stripped):
            block_lines = [line]
            tag_m = re.match(r"^<(p|h[1-6]|div|details)\b", stripped)
            if tag_m and f"</{tag_m.group(1)}>" not in stripped:
                end_tag = f"</{tag_m.group(1)}>"
                i += 1
                while i < n:
                    block_lines.append(lines[i])
                    if end_tag in lines[i]:
                        break
                    i += 1
            raw_block = "\n".join(block_lines)
            raw_block = re.sub(
                r'(\b(?:href|src)=")([^"]+)(")',
                lambda am: am.group(1) + resolve_href(html.unescape(am.group(2)), source_md, repo_root) + am.group(3),
                raw_block,
            )
            out.append(raw_block)
            i += 1
            continue

        # Fenced code block
        fence_m = re.match(r"^```(\w*)\s*$", stripped)
        if fence_m:
            lang = fence_m.group(1)
            code_lines: list[str] = []
            i += 1
            while i < n and not lines[i].strip().startswith("```"):
                code_lines.append(lines[i])
                i += 1
            cls = f' class="language-{html.escape(lang, quote=True)}"' if lang else ""
            out.append(f"<pre><code{cls}>{html.escape(chr(10).join(code_lines), quote=False)}</code></pre>")
            i += 1
            continue

        # Horizontal rule
        if re.fullmatch(r"-{3,}|\*{3,}", stripped):
            out.append("<hr />")
            i += 1
            continue

        # ATX Heading
        h_m = re.match(r"^(#{1,6})\s+(.+?)\s*$", line)
        if h_m:
            level = len(h_m.group(1))
            raw_title = h_m.group(2)
            base_slug = github_slug(raw_title)
            count = slug_counts.get(base_slug, 0)
            slug_counts[base_slug] = count + 1
            slug = base_slug if count == 0 else f"{base_slug}-{count}"
            rendered_title = render_inline(raw_title, source_md, repo_root)
            out.append(f'<h{level} id="{html.escape(slug, quote=True)}">{rendered_title}</h{level}>')
            if level in (2, 3):
                plain = html.unescape(re.sub(r"<[^>]+>", "", rendered_title))
                toc.append((level, slug, plain))
            i += 1
            continue

        # Pipe Table
        if "|" in line and i + 1 < n and is_table_divider(lines[i + 1]):
            headers = split_table_row(line)
            aligns_raw = split_table_row(lines[i + 1])
            aligns: list[str] = []
            for a in aligns_raw:
                if a.startswith(":") and a.endswith(":"):
                    aligns.append(' style="text-align:center"')
                elif a.endswith(":"):
                    aligns.append(' style="text-align:right"')
                else:
                    aligns.append("")
            i += 2
            rows: list[list[str]] = []
            while i < n and "|" in lines[i] and lines[i].strip():
                rows.append(split_table_row(lines[i]))
                i += 1
            th_html = "".join(
                f"<th{aligns[idx] if idx < len(aligns) else ''}>{render_inline(cell, source_md, repo_root)}</th>"
                for idx, cell in enumerate(headers)
            )
            body_rows: list[str] = []
            for r in rows:
                tds = "".join(
                    f"<td{aligns[idx] if idx < len(aligns) else ''}>{render_inline(cell, source_md, repo_root)}</td>"
                    for idx, cell in enumerate(r)
                )
                body_rows.append(f"<tr>{tds}</tr>")
            out.append(f"<table><thead><tr>{th_html}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>")
            continue

        # Blockquote / GitHub Alert
        if stripped.startswith(">"):
            q_lines: list[str] = []
            while i < n and lines[i].strip().startswith(">"):
                q_lines.append(re.sub(r"^\s*>\s?", "", lines[i]))
                i += 1
            alert_m = re.match(r"^\[!(NOTE|TIP|IMPORTANT|WARNING|CAUTION)\]\s*$", q_lines[0].strip(), re.I)
            if alert_m:
                kind = alert_m.group(1).lower()
                body_html, _ = render_markdown("\n".join(q_lines[1:]), source_md, repo_root)
                out.append(
                    f'<div class="gh-alert gh-alert-{kind}">'
                    f'<div class="gh-alert-title">{kind.upper()}</div>'
                    f"{body_html}</div>"
                )
            else:
                body_html, _ = render_markdown("\n".join(q_lines), source_md, repo_root)
                out.append(f"<blockquote>{body_html}</blockquote>")
            continue

        # Lists (ordered or unordered, supporting indented sub-bullets)
        list_m = re.match(r"^(\s*)(?:([-*])|(\d+)\.)\s+(.+)$", line)
        if list_m:
            ordered = bool(list_m.group(3))
            tag = "ol" if ordered else "ul"
            items_html: list[str] = []
            while i < n:
                m_item = re.match(r"^(\s*)(?:([-*])|(\d+)\.)\s+(.+)$", lines[i])
                if not m_item or len(m_item.group(1)) > 1 or (bool(m_item.group(3)) != ordered):
                    break
                main_text = render_inline(m_item.group(4), source_md, repo_root)
                sub_items: list[str] = []
                i += 1
                while i < n:
                    m_sub = re.match(r"^(\s{2,})(?:[-*]|\d+\.)\s+(.+)$", lines[i])
                    if not m_sub:
                        break
                    sub_items.append(f"<li>{render_inline(m_sub.group(2), source_md, repo_root)}</li>")
                    i += 1
                sub_html = f"<ul>{''.join(sub_items)}</ul>" if sub_items else ""
                items_html.append(f"<li>{main_text}{sub_html}</li>")
            out.append(f"<{tag}>{''.join(items_html)}</{tag}>")
            continue

        # Paragraph
        p_lines = [stripped]
        i += 1
        while i < n:
            nxt = lines[i].strip()
            if (
                not nxt
                or nxt.startswith(("#", ">", "```", "<p", "<h", "<div", "<a id="))
                or re.fullmatch(r"-{3,}|\*{3,}", nxt)
                or re.match(r"^\s*(?:[-*]|\d+\.)\s+", lines[i])
                or ("|" in lines[i] and i + 1 < n and is_table_divider(lines[i + 1]))
            ):
                break
            p_lines.append(nxt)
            i += 1
        out.append(f"<p>{render_inline(' '.join(p_lines), source_md, repo_root)}</p>")

    return "\n".join(out), toc


def wrap_doc_page(title: str, lang: str, active_href: str, body_html: str, toc: list[tuple[int, str, str]]) -> str:
    toc_items = "".join(
        f'<li style="padding-left:{0.75 if lvl == 3 else 0}rem"><a href="#{html.escape(slug, quote=True)}">{html.escape(label)}</a></li>'
        for lvl, slug, label in toc
    )
    toc_heading = "目次" if lang == "ja" else "On this page"
    toc_aside = f'<aside class="doc-toc"><h4>{toc_heading}</h4><ul>{toc_items}</ul></aside>' if toc_items else ""

    nav_items = [
        ("index.html", "Home"),
        ("desktop-guide.html", "Desktop Guide (EN)"),
        ("desktop-guide-ja.html", "日本語ガイド"),
        ("analysis-scope.html", "Interpretation Scope"),
        ("manuscript-guide.html", "Manuscript Guide"),
    ]
    nav_links = "".join(
        f'<a href="{href}"{" class=\"active\"" if href == active_href else ""}>{label}</a>'
        for href, label in nav_items
    )

    return f"""<!DOCTYPE html>
<html lang="{lang}">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(title)}</title>
  <link rel="icon" type="image/svg+xml" href="assets/images/mfrm-desktop-icon.svg" />
  <link rel="stylesheet" href="styles.css" />
</head>
<body>
  <header class="site-nav">
    <div class="container nav-inner">
      <a class="brand-link" href="index.html">
        <img src="assets/images/mfrm-desktop-icon.svg" alt="MFRM icon" />
        <span>MFRM</span>
      </a>
      <nav class="nav-links" aria-label="Documentation">
        {nav_links}
      </nav>
      <div class="nav-actions">
        <a class="nav-gh" href="https://mfrmpythonapplication-xodghmcghngx8cz5zmhfeg.streamlit.app/" target="_blank" rel="noopener">Web App ↗</a>
        <a class="nav-gh" href="https://ryuya-dot-com.github.io/mfrmr/" target="_blank" rel="noopener">mfrmr (R) ↗</a>
        <a class="nav-gh" href="https://github.com/Ryuya-dot-com/MFRM_Python_Application" target="_blank" rel="noopener">GitHub ↗</a>
      </div>
    </div>
  </header>

  <main class="container doc-page-wrap">
    <article class="doc-article">
      {body_html}
    </article>
    {toc_aside}
  </main>

  <footer class="site-footer">
    <div class="container footer-inner">
      <div>
        <strong>MFRM</strong> · <a href="index.html">Product Overview</a> ·
        <a href="desktop-guide.html">Desktop Guide</a> ·
        <a href="analysis-scope.html">Interpretation Scope</a>
      </div>
      <div>
        <a href="https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-windows-20261030">Windows Preview</a> ·
        <a href="https://mfrmpythonapplication-xodghmcghngx8cz5zmhfeg.streamlit.app/" target="_blank" rel="noopener">Streamlit Web</a> ·
        <a href="https://ryuya-dot-com.github.io/mfrmr/" target="_blank" rel="noopener">mfrmr (R)</a>
      </div>
    </div>
  </footer>
</body>
</html>
"""


def verify_output(output_dir: Path, repo_root: Path) -> None:
    html_files = sorted(output_dir.glob("*.html"))
    assert html_files, "No HTML files generated"

    # Collect all id="..." anchors per HTML file
    anchors_by_file: dict[Path, set[str]] = {}
    for hf in html_files:
        txt = hf.read_text(encoding="utf-8")
        anchors_by_file[hf.resolve()] = set(re.findall(r'\bid="([^"]+)"', txt))

    # Verify all local href/src targets & internal #anchors, plus zero unrendered markdown links
    attr_re = re.compile(r'\b(?:href|src)="([^"]+)"')
    workstation_path_re = re.compile(r"(?:/Users/|/private/tmp|C:\\Users\\)")

    for hf in html_files:
        txt = hf.read_text(encoding="utf-8")
        m_path = workstation_path_re.search(txt)
        assert m_path is None, f"Workstation path {m_path.group(0)!r} in {hf.name}"
        for raw in attr_re.findall(txt):
            val = html.unescape(raw).strip()
            if not val or re.match(r"^[a-zA-Z][a-zA-Z0-9+.-]*:|^//|^mailto:", val):
                continue
            path_part, sep, frag = val.partition("#")
            target = (hf.parent / path_part).resolve() if path_part else hf.resolve()
            assert target.exists(), f"Broken local reference {val!r} in {hf.name}"
            if sep and target.suffix == ".html":
                avail = anchors_by_file.get(target, set())
                assert frag in avail, f"Missing anchor #{frag} in {target.name} (from {hf.name})"

    # Respectful-language check on landing page source
    respectful_patterns = [
        re.compile(r"\bbeginners?(?:[- ](?:friendly|oriented|facing))?\b", re.I),
        re.compile(r"\bnovices?\b", re.I),
        re.compile(r"\bnon[- ]?experts?\b", re.I),
        re.compile(r"\bnon[- ]?technical(?:\s+(?:users?|audience|readers?))?\b", re.I),
        re.compile(r"\blay(?:persons?|people|users?|audience|readers?)\b", re.I),
        re.compile(r"\binexperienced\b", re.I),
        re.compile(r"\bless experienced\b", re.I),
        re.compile(r"初心者|初学者|ビギナー|素人|非専門家"),
    ]
    landing_txt = (repo_root / "site" / "index.html").read_text(encoding="utf-8")
    for pat in respectful_patterns:
        m = pat.search(landing_txt)
        assert m is None, f"Forbidden ability label {m.group(0)!r} in site/index.html"


def build_site(repo_root: Path, output_dir: Path) -> None:
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / ".nojekyll").write_text("", encoding="utf-8")

    # Copy static landing assets
    shutil.copy2(repo_root / "site" / "index.html", output_dir / "index.html")
    shutil.copy2(repo_root / "site" / "styles.css", output_dir / "styles.css")

    img_dst = output_dir / "assets" / "images"
    img_dst.mkdir(parents=True, exist_ok=True)
    for img in sorted((repo_root / "docs" / "images").iterdir()):
        if img.is_file():
            shutil.copy2(img, img_dst / img.name)

    ex_dst = output_dir / "examples"
    ex_dst.mkdir(parents=True, exist_ok=True)
    shutil.copy2(repo_root / "examples" / "desktop_writing.csv", ex_dst / "desktop_writing.csv")

    # Render canonical Markdown guides
    for rel_md, out_name, lang, page_title in PAGE_SPECS:
        src_path = repo_root / rel_md
        body_html, toc = render_markdown(src_path.read_text(encoding="utf-8"), src_path, repo_root)
        page_html = wrap_doc_page(page_title, lang, out_name, body_html, toc)
        (output_dir / out_name).write_text(page_html, encoding="utf-8")

    # Clean 404 page
    not_found = wrap_doc_page(
        "Page Not Found — MFRM",
        "en",
        "",
        '<h1>Page Not Found</h1><p>Return to the <a href="index.html">MFRM homepage</a> or open the <a href="desktop-guide.html">Desktop User Guide</a>.</p>',
        [],
    )
    (output_dir / "404.html").write_text(not_found, encoding="utf-8")

    verify_output(output_dir, repo_root)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="Directory to write the built static site")
    args = parser.parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    build_site(repo_root, args.output.resolve())
    print(f"Built and verified GitHub Pages site at {args.output.resolve()}")


if __name__ == "__main__":
    main()

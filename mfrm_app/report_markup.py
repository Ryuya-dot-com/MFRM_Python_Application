"""The Markdown subset emitted by the methods and writing appendices."""

import html
import re


def inline_spans(text):
    """Keep literal text safe while recognizing emphasis and inline code."""
    pattern = r"\\([\\`*|])|`([^`\n]+)`|\*\*(.+?)\*\*|(?<!\w)\*([^*\n]+)\*(?!\w)"
    end = 0
    for match in re.finditer(pattern, text):
        if match.start() > end:
            yield "text", text[end:match.start()]
        for kind, value in zip(("text", "code", "bold", "italic"), match.groups()):
            if value is not None:
                yield kind, value
                break
        end = match.end()
    if end < len(text):
        yield "text", text[end:]


def inline_html(text):
    tags = {"bold": "b", "italic": "i", "code": "code"}
    return "".join(
        f"<{tags[kind]}>{html.escape(value)}</{tags[kind]}>" if kind in tags else html.escape(value)
        for kind, value in inline_spans(text)
    )


def blocks(markdown):
    """Yield (kind, depth, content); preserve unknown syntax as literal text.

    shortcut: supports the generated headings, bullets and pipe tables, not full
    CommonMark; extend against real appendix examples if generators add syntax.
    """
    lines = markdown.splitlines()

    def cells(line):
        return [cell.strip() for cell in re.split(r"(?<!\\)\|", line.strip()[1:-1])]

    i = 0
    while i < len(lines):
        raw = lines[i]
        text = raw.strip()
        i += 1
        if not text:
            continue
        if text == "---":
            yield "rule", 0, ""
        elif heading := re.match(r"^(#{1,3})\s+(.+)$", text):
            yield "heading", len(heading[1]), heading[2]
        elif text.startswith(("- ", "* ")):
            yield "bullet", (len(raw) - len(raw.lstrip())) // 2, text[2:]
        elif (text.startswith("|") and text.endswith("|") and i < len(lines)
              and lines[i].strip().startswith("|") and lines[i].strip().endswith("|")
              and all(re.fullmatch(r":?-{3,}:?", cell) for cell in cells(lines[i]))):
            rows = [cells(text)]
            i += 1  # Only the alignment row is markup; keep blank and dash-only data cells.
            while i < len(lines) and lines[i].strip().startswith("|") and lines[i].strip().endswith("|"):
                rows.append(cells(lines[i]))
                i += 1
            width = max(map(len, rows))
            yield "table", 0, [row + [""] * (width - len(row)) for row in rows]
        else:
            yield "paragraph", 0, text

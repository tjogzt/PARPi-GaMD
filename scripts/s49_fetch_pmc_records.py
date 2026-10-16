#!/usr/bin/env python3
#
# Purpose:  extract quantitative trapping data from PMC JATS full-text XML
# Created:  2026-09 (header standardised 2026-10-11)
# Run:      python3 scripts/s49_fetch_pmc_records.py   (from the repository root)

"""
extract_pmc.py -- extract quantitative trapping data from PMC JATS full-text XML

Extracts: title/abstract, all tables (as text), and sentences containing trapping/retention keywords with numeric values.
Outputs markdown for manual verification and inclusion in the trapping dataset (quantitative fold-trapping / EC50 etc.).

Usage: python extract_pmc.py data/00_raw/literature/PMC3528345.xml
"""
from __future__ import annotations
import re
import sys
from pathlib import Path
import xml.etree.ElementTree as ET

KEYWORDS = re.compile(
    r"\b(trap|trapp|retention|retain|PARP[- ]?DNA|chromatin[- ]?bound|"
    r"allosteric|fold|EC50|IC50|Kd|dissociat|cytotox)\b", re.I)
HAS_NUM = re.compile(r"\d")


def txt(el) -> str:
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip()


def split_sentences(s: str) -> list[str]:
    return re.split(r"(?<=[.;])\s+(?=[A-Z(])", s)


def main() -> int:
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    xml = Path(sys.argv[1])
    root = ET.fromstring(xml.read_text(errors="ignore"))
    out = []

    title = root.find(".//article-title")
    out.append(f"# PMC extract: {txt(title) if title is not None else xml.stem}\n")

    abs = root.find(".//abstract")
    if abs is not None:
        out.append("## Abstract\n" + txt(abs) + "\n")

    # tables
    tables = root.findall(".//table-wrap")
    out.append(f"\n## Tables ({len(tables)})\n")
    for i, tw in enumerate(tables, 1):
        cap = tw.find(".//caption")
        label = tw.find(".//label")
        out.append(f"\n### Table {i}: "
                   f"{txt(label) if label is not None else ''} "
                   f"{txt(cap) if cap is not None else ''}")
        for tr in tw.findall(".//tr"):
            cells = [txt(c) for c in tr.findall("./*")]
            if any(cells):
                out.append("| " + " | ".join(cells) + " |")

    # body-text key sentences
    out.append("\n## Sentences with trapping keywords + numeric values\n")
    seen = set()
    for p in root.findall(".//body//p"):
        for sent in split_sentences(txt(p)):
            if KEYWORDS.search(sent) and HAS_NUM.search(sent) and sent not in seen:
                seen.add(sent)
                out.append(f"- {sent}")

    # figure captions
    out.append("\n## Figure captions containing keywords\n")
    for fig in root.findall(".//fig"):
        c = fig.find(".//caption")
        if c is not None:
            t = txt(c)
            if KEYWORDS.search(t):
                lab = fig.find(".//label")
                out.append(f"- **{txt(lab) if lab is not None else 'Fig'}**: {t}")

    dest = xml.with_name(xml.stem + "_extracted.md")
    dest.write_text("\n".join(out) + "\n")
    print(f"[OK] -> {dest}  ({len(tables)} tables, {len(seen)} key sentences)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

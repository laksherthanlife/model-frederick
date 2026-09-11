"""Repair a Bio-Rad export whose zip entries were written in lowercase.

The 2026-07-24 qPCR export cannot be opened by any standard reader: its archive
names are lowercased (`[content_types].xml`, `xl/sharedstrings.xml`), and the OOXML
readers match those case-sensitively. Nothing is wrong with the data.

Usage: python scripts/repair_qpcr_export.py <in.xlsx> <out.xlsx>
"""
from __future__ import annotations

import pathlib
import re
import sys
import zipfile

CANONICAL = {
    "[content_types].xml": "[Content_Types].xml",
    "_rels/.rels": "_rels/.rels",
    "xl/workbook.xml": "xl/workbook.xml",
    "xl/sharedstrings.xml": "xl/sharedStrings.xml",
    "xl/styles.xml": "xl/styles.xml",
    "xl/_rels/workbook.xml.rels": "xl/_rels/workbook.xml.rels",
}


def canonical_name(name: str) -> str:
    normalised = name.replace("\\", "/")
    lowered = normalised.lower()
    if lowered in CANONICAL:
        return CANONICAL[lowered]
    sheet = re.fullmatch(r"xl/worksheets/sheet(\d+)\.xml", lowered)
    return f"xl/worksheets/sheet{sheet.group(1)}.xml" if sheet else normalised


def repair(src: pathlib.Path, dst: pathlib.Path) -> pathlib.Path:
    dst.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(src) as zin, zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            name = canonical_name(item.filename)
            # Relationship targets inside the XML are lowercased too.
            if name.endswith((".rels", "workbook.xml")):
                data = data.replace(b"sharedstrings.xml", b"sharedStrings.xml")
            zout.writestr(name, data)
    return dst


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    print("repaired ->", repair(pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])))

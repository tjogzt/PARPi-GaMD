#!/usr/bin/env python3
"""verify_manifest.py — cross-check data_manifest.md checksums against files on disk.

Parses the manifest table (data_manifest.md), and for every row whose status is
'current' and whose checksum is a real md5 (32 hex chars), verifies the artifact's
current md5. Exits non-zero on any mismatch. Rows marked 'legacy', 'superseded',
or '(see notes)' are reported but not enforced.

Usage: python3 scripts/verify_manifest.py [--root <repo root>]
"""
import argparse
import hashlib
import re
import sys
from pathlib import Path

MD5_RE = re.compile(r"^[0-9a-f]{32}$")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    args = ap.parse_args()
    root = Path(args.root).resolve()

    manifest = root / "data_manifest.md"
    if not manifest.exists():
        print(f"manifest not found: {manifest}")
        return 2

    rows = []
    for line in manifest.read_text().splitlines():
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 10 or cells[0] == "artifact":
            continue
        artifact, producer, line_no, upstream, downstream, ms_loc, protocol, date, status, checksum = cells[:10]
        rows.append((artifact, status, checksum))

    ok, checked, skipped, missing, mismatch = 0, 0, 0, [], []
    for artifact, status, checksum in rows:
        if status != "current":
            skipped += 1
            continue
        if not MD5_RE.match(checksum):
            skipped += 1  # '(see notes)' etc.
            continue
        checked += 1
        f = root / artifact
        if not f.exists():
            missing.append(artifact)
            continue
        actual = hashlib.md5(f.read_bytes()).hexdigest()
        if actual == checksum:
            ok += 1
        else:
            mismatch.append((artifact, checksum, actual))

    print(f"verified {checked} checksums: {ok} OK, {len(mismatch)} MISMATCH, {len(missing)} MISSING, {skipped} skipped")
    for a, exp, act in mismatch:
        print(f"  MISMATCH {a}\n    expected {exp}\n    actual   {act}")
    for a in missing:
        print(f"  MISSING  {a}")
    return 1 if (mismatch or missing) else 0


if __name__ == "__main__":
    sys.exit(main())

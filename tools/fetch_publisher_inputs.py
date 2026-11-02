#!/usr/bin/env python3
"""Retrieve publisher inputs only when they match the recorded SHA-256 digests.

Uses only the Python standard library. Existing matching files are retained.
Existing mismatched files are never overwritten. No reference edition is substituted.
"""
from __future__ import annotations
import argparse
import hashlib
from pathlib import Path
import sys
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
BASE = "https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41567-026-03442-w/MediaObjects/"

def retrieve(relative: str, url: str, expected: str) -> None:
    target = ROOT / relative
    if target.is_file():
        actual = hashlib.sha256(target.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"Existing file differs from recorded input: {relative}")
        print(f"Verified existing: {relative}")
        return
    request = urllib.request.Request(url, headers={"User-Agent": "Thermodynamic-Accessibility-reproduction/4"})
    with urllib.request.urlopen(request, timeout=90) as response:
        data = response.read()
    actual = hashlib.sha256(data).hexdigest()
    if actual != expected:
        raise ValueError(f"Downloaded edition differs: {relative}; expected {expected}, received {actual}")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".part")
    temporary.write_bytes(data)
    temporary.replace(target)
    print(f"Downloaded and verified: {relative}")

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--with-papers", action="store_true", help="Also retrieve the three publisher PDFs when their exact bytes are available.")
    args = parser.parse_args()
    entries = {}
    for line in (ROOT / "REFERENCE_MANIFEST.sha256").read_text(encoding="utf-8").splitlines():
        digest, relative = line.split(maxsplit=1)
        entries[relative] = digest
    failures = []
    for number in range(3, 7):
        name = f"41567_2026_3442_MOESM{number}_ESM.xlsx"
        relative = "data/" + name
        try:
            retrieve(relative, BASE + name, entries[relative])
        except Exception as exc:
            failures.append(relative)
            print(f"NOT RETRIEVED: {relative}: {exc}", file=sys.stderr)
    if args.with_papers:
        sources = {
            "references/s41567-026-03442-w.pdf": "https://www.nature.com/articles/s41567-026-03442-w.pdf",
            "references/41567_2026_3442_MOESM1_ESM.pdf": BASE + "41567_2026_3442_MOESM1_ESM.pdf",
            "references/41567_2026_3442_MOESM2_ESM.pdf": BASE + "41567_2026_3442_MOESM2_ESM.pdf",
        }
        for relative, url in sources.items():
            try:
                retrieve(relative, url, entries[relative])
            except Exception as exc:
                failures.append(relative)
                print(f"NOT RETRIEVED: {relative}: {exc}", file=sys.stderr)
    if failures:
        print("Use the exact files from the companion distribution for unresolved entries; no mismatched files were accepted.", file=sys.stderr)
        return 1
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Fetch verified MOR comparator structures from the RCSB PDB into data/comparators/.

Comparator set and rationale are documented in docs/references.md. Idempotent: skips
files already present. Run: ``python src/01.01.00_fetch_comparators.py``.
"""

from __future__ import annotations

import sys
import urllib.error
import urllib.request

# Add src/ to path so the numbered script can import the shared package.
sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

from zh853mor import paths  # noqa: E402

RCSB_URL = "https://files.rcsb.org/download/{pdb}.pdb"

# (PDB id, short description) -- see docs/references.md for full citations.
COMPARATORS: list[tuple[str, str]] = [
    ("8F7R", "endomorphin-1 / Gi / scFv16 (closest chemical analog)"),
    ("8EFQ", "DAMGO / Gi (human)"),
    ("6DDE", "DAMGO / Gi (mouse, 3.5 A)"),
    ("8F7Q", "beta-endorphin / Gi"),
    ("5C1M", "BU72 / Nb39 (2.07 A X-ray, high-res pocket)"),
    ("7T2G", "mitragynine pseudoindoxyl / Gi (2.5 A, biased)"),
    ("8EF5", "fentanyl / Gi"),
    ("8EFB", "oliceridine (TRV130) / Gi (biased)"),
    ("8EFL", "SR-17018 / Gi (biased)"),
    ("8EFO", "PZM21 / Gi (biased)"),
    ("4DKL", "beta-FNA (inactive/antagonist reference)"),
    ("9WST", "DAMGO / Gz (bias comparator)"),
    ("9WSV", "DAMGO / beta-arrestin-1 (bias comparator)"),
]


def fetch(pdb: str, dest_dir: str) -> bool:
    """Download ``pdb`` into ``dest_dir``; return True if fetched, False if skipped."""
    out = paths.COMPARATORS / f"{pdb}.pdb"
    if out.exists():
        print(f"  skip {pdb} (exists)")
        return False
    url = RCSB_URL.format(pdb=pdb)
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:  # noqa: S310 (trusted RCSB host)
            out.write_bytes(resp.read())
    except (urllib.error.URLError, TimeoutError) as exc:
        print(f"  FAIL {pdb}: {exc}", file=sys.stderr)
        return False
    print(f"  got  {pdb}  ({out.stat().st_size // 1024} KB)  {dict(COMPARATORS)[pdb]}")
    return True


def main() -> int:
    paths.ensure_dir(paths.COMPARATORS)
    print(f"Fetching {len(COMPARATORS)} comparator structures into {paths.COMPARATORS}")
    fetched = sum(fetch(pdb, str(paths.COMPARATORS)) for pdb, _ in COMPARATORS)
    print(f"Done: {fetched} fetched, {len(COMPARATORS) - fetched} skipped/failed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

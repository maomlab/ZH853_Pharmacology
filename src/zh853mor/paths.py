"""Canonical project paths.

Resolved relative to the repository root so scripts work regardless of the current
working directory. Follows the OBJECTIVES.md layout: data/ src/ intermediate/ product/.
"""

from __future__ import annotations

from pathlib import Path

# src/zh853mor/paths.py -> repo root is three parents up.
ROOT: Path = Path(__file__).resolve().parents[2]

DATA: Path = ROOT / "data"
SRC: Path = ROOT / "src"
INTERMEDIATE: Path = ROOT / "intermediate"
PRODUCT: Path = ROOT / "product"
DOCS: Path = ROOT / "docs"

# The deposited cryo-EM model (MOR-Gi-scFv16-ZH853).
CRYOEM_PDB: Path = DATA / "mor_gi_scfv16_refine-coot-22_real_space_refined_169_edit.pdb"

# Downloaded comparator structures land here.
COMPARATORS: Path = DATA / "comparators"


def ensure_dir(path: Path) -> Path:
    """Create ``path`` (and parents) if absent and return it."""
    path.mkdir(parents=True, exist_ok=True)
    return path

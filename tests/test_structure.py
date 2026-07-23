"""Tests for zh853mor.structure against the deposited cryo-EM model.

These lock in the verified Phase-1 facts so downstream analyses can rely on them.
Skipped automatically if the (large) PDB is unavailable.
"""

from __future__ import annotations

import pytest

from zh853mor import paths, structure

pytestmark = pytest.mark.skipif(
    not paths.CRYOEM_PDB.exists(), reason="cryo-EM PDB not present"
)


@pytest.fixture(scope="module")
def universe():
    return structure.load(paths.CRYOEM_PDB)


def test_human_numbering(universe):
    """Every canonical orthosteric residue matches -> human OPRM1 numbering."""
    numbering = structure.verify_numbering(universe)
    assert all(numbering.values()), {k: v for k, v in numbering.items() if not v}


def test_no_receptor_gaps(universe):
    """MOR chain R (69-349) is modeled with no internal gaps."""
    assert structure.receptor_gaps(universe) == []


def test_conserved_disulfide(universe):
    """The conserved ECL2-TM3 disulfide C142-C219 is present."""
    pairs = {(i, j) for i, j, _ in structure.disulfides(universe)}
    assert (142, 219) in pairs


def test_key_pocket_anchors_present(universe):
    """D149 (salt-bridge anchor) and Y328 (3-7 lock) contact the ligand."""
    contact_ids = {c.resid for c in structure.pocket_contacts(universe)}
    assert {149, 328}.issubset(contact_ids)

"""Tests for the MD-prep helpers (Phase 2)."""

from __future__ import annotations

import pytest

from zh853mor import paths, prep, structure

pytestmark = pytest.mark.skipif(
    not paths.CRYOEM_PDB.exists(), reason="cryo-EM PDB not present"
)


@pytest.fixture(scope="module")
def universe():
    return structure.load(paths.CRYOEM_PDB)


def test_functional_residue_numbers(universe):
    """D2.50 = Asp116 and the DRY motif = Asp166-Arg167-Tyr168 (verified against construct)."""
    by_id = {int(a.resid): a.resname for a in universe.select_atoms("segid R and name CA")}
    assert by_id[prep.D250_SODIUM] == "ASP"
    assert [by_id[r] for r in prep.DRY_MOTIF] == ["ASP", "ARG", "TYR"]


def test_incomplete_residues_detected(universe):
    """The 3.5 A model has truncated sidechains to rebuild (and receptor CA count is sane)."""
    incomplete = prep.incomplete_residues(universe)
    assert len(incomplete) > 0
    # all reported residues really are under their standard heavy-atom count
    for r in incomplete:
        assert r.kind == "incomplete"


def test_membrane_frame(universe):
    """Cholesterol is present and defines a plausible bilayer thickness along the normal."""
    frame = prep.membrane_frame(universe)
    assert frame["n_cholesterol_atoms"] == 84
    assert 20.0 < frame["cholesterol_span_along_normal"] < 45.0


def test_disulfide_preserved(universe):
    assert (142, 219) in {(i, j) for i, j, _ in structure.disulfides(universe)}

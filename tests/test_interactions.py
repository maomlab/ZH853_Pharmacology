"""Tests for the interaction-fingerprint engine and comparator loading.

Lock in the key Objective-1/2 findings so refactors can't silently change them.
Skipped if the comparator PDBs have not been fetched.
"""

from __future__ import annotations

import pytest

from zh853mor import comparators, interactions, paths

pytestmark = pytest.mark.skipif(
    not (paths.COMPARATORS / "8F7R.pdb").exists() or not paths.CRYOEM_PDB.exists(),
    reason="comparator PDBs not fetched",
)


@pytest.fixture(scope="module")
def zh853_fp():
    cx = comparators.load_complex("ZH853")
    return interactions.fingerprint(cx.receptor, cx.ligand, cx.offset_to_human)


def test_offset_detection():
    """Human structures map with offset 0; mouse structures with +2."""
    assert comparators.load_complex("8EFQ").offset_to_human == 0  # human DAMGO
    assert comparators.load_complex("5C1M").offset_to_human == 2  # mouse BU72


def test_zh853_salt_bridges(zh853_fp):
    """ZH853 makes ionic contacts at the D3.32 anchor and the distinctive ECL2 E231."""
    assert "ionic" in zh853_fp[149].interactions
    assert "ionic" in zh853_fp[231].interactions


def test_zh853_aromatic_ring_detection(zh853_fp):
    """Geometric ring perception finds the H321 aromatic/cation-pi contact."""
    assert zh853_fp[321].interactions & {"aromatic", "cation_pi"}


def test_e231_is_zh853_distinctive():
    """No other agonist contacts ECL2 Glu231 — the key selectivity lead."""
    others = ["8F7R", "8EFQ", "6DDE", "8F7Q", "5C1M", "8EF5", "8EFB", "8EFL", "8EFO", "7T2G"]
    for pdb in others:
        cx = comparators.load_complex(pdb)
        fp = interactions.fingerprint(cx.receptor, cx.ligand, cx.offset_to_human)
        assert 231 not in fp, f"{pdb} unexpectedly contacts residue 231"


def test_universal_anchor_shared():
    """D149 (D3.32) is contacted by every agonist (universal anchor, not selective)."""
    for pdb in ["8F7R", "8EFQ", "5C1M", "8EF5", "8EFO"]:
        cx = comparators.load_complex(pdb)
        fp = interactions.fingerprint(cx.receptor, cx.ligand, cx.offset_to_human)
        assert 149 in fp

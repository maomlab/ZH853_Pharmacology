"""Tests for the cheminformatics module (Objective 3)."""

from __future__ import annotations

from zh853mor import chem


def test_all_analogs_parse():
    for name in chem.ANALOGS:
        assert chem.mol(name).GetNumHeavyAtoms() > 40


def test_zh853_matches_deposited_ligand():
    """ZH853 SMILES has 59 heavy atoms — matches the deposited L01 ligand."""
    assert chem.mol("ZH853").GetNumHeavyAtoms() == 59


def test_all_analogs_are_bro5():
    """Every analog violates rule-of-5 (the core Objective-3 liability)."""
    for name in chem.ANALOGS:
        p = chem.descriptors_2d(name, chem.mol(name))
        assert "HBD>5" in p.bro5_flags and "TPSA>140" in p.bro5_flags


def test_n_methylation_reduces_donors():
    """Each N-methyl removes one H-bond donor and lowers TPSA."""
    m = chem.mol("ZH853")
    p0 = chem.descriptors_2d("ZH853", m)
    sites = chem.n_methyl_sites(m)
    assert len(sites) == 6
    tri = chem.add_n_methyls(m, sites[:3])
    p = chem.descriptors_2d("tri", tri)
    assert p.hbd == p0.hbd - 3
    assert p.tpsa < p0.tpsa


def test_lipidation_builds_valid_mol():
    """Semaglutide-style conjugation yields a valid, heavier, more complex molecule."""
    m = chem.mol("ZH853")
    conj = chem.conjugate_at_primary_amide(m, chem.LINKERS["palmitoyl"][0])
    p0 = chem.descriptors_2d("ZH853", m)
    p = chem.descriptors_2d("palmitoyl", conj)
    assert p.mw > p0.mw + 150
    assert p.clogp > p0.clogp  # lipidation raises lipophilicity

"""Tests for zh853mor.md on a synthetic membrane system.

No trajectory from the cluster is needed (or wanted) here: the things that can silently go wrong
-- transferring the construct numbering onto a tleap-renumbered topology, reducing a distance
matrix per residue, counting a Lipid21 POPC as one lipid rather than three -- are all exercised
by a small system whose right answers are known by construction.
"""

from __future__ import annotations

import MDAnalysis as mda
import numpy as np
import pytest

from zh853mor import md

# A miniature receptor: the construct numbers are 69.. as in the real staged receptor.pdb, while
# the built system is renumbered from 1 by tleap -- which is exactly the mismatch md.py exists to
# resolve.
RESNAMES = ["ALA", "ASP", "GLU", "TYR", "HIS", "ASN"]
FIRST_CONSTRUCT_RESID = 69


def _universe(resnames, resids, names, masses, positions, segid="SYSTEM"):
    n = len(names)
    resindex = np.repeat(np.arange(len(resnames)), [sum(1 for r in resids if r == x)
                                                    for x in sorted(set(resids), key=resids.index)])
    u = mda.Universe.empty(n, n_residues=len(resnames), atom_resindex=resindex,
                           residue_segindex=np.zeros(len(resnames), dtype=int), trajectory=True)
    u.add_TopologyAttr("name", names)
    u.add_TopologyAttr("resname", resnames)
    u.add_TopologyAttr("resid", sorted(set(resids), key=resids.index))
    u.add_TopologyAttr("segid", [segid])
    u.add_TopologyAttr("mass", masses)
    u.load_new(np.asarray(positions, dtype=np.float32)[None, :, :], order="fac")
    return u


def _protein_atoms(resnames, start_resid, x0=0.0):
    """Backbone N/CA/C/O per residue, 4 A apart along x."""
    names, resids, masses, pos = [], [], [], []
    for i, _ in enumerate(resnames):
        for name, mass in (("N", 14.01), ("CA", 12.01), ("C", 12.01), ("O", 16.00)):
            names.append(name)
            resids.append(start_resid + i)
            masses.append(mass)
            pos.append([x0 + 4.0 * i, 0.0, 0.0])
    return names, resids, masses, pos


@pytest.fixture
def receptor_pdb(tmp_path):
    """The staged receptor, in construct numbering (69...)."""
    names, resids, masses, pos = _protein_atoms(RESNAMES, FIRST_CONSTRUCT_RESID)
    u = _universe(RESNAMES, resids, names, masses, pos, segid="R")
    path = tmp_path / "receptor.pdb"
    u.atoms.write(str(path))
    return path


@pytest.fixture
def system():
    """The built system: protein renumbered from 1, plus ligand, water, Na+ and lipids."""
    names, resids, masses, pos = _protein_atoms(RESNAMES, 1)
    resnames = list(RESNAMES)

    def add(resname, atoms):
        """Append one residue; `atoms` is [(name, mass, position), ...]."""
        resnames.append(resname)
        for name, mass, xyz in atoms:
            names.append(name)
            resids.append(len(resnames))
            masses.append(mass)
            pos.append(xyz)

    # Ligand: two polar atoms, 3.0 A from residue 3 (construct 71) and further from every other.
    add("LIG", [("N01", 14.01, [8.0, 3.0, 0.0]), ("O02", 16.00, [8.0, 3.4, 0.0])])
    # One bridging water: within 3.5 A of the ligand and of residue 3 (construct 71) only.
    add("WAT", [("O", 16.00, [8.0, 2.5, 1.5]), ("H1", 1.008, [8.5, 2.5, 2.0]),
                ("H2", 1.008, [7.5, 2.5, 2.0])])
    add("Na+", [("Na+", 22.99, [50.0, 50.0, 0.0])])

    # Four POPC (Lipid21: PC + PA + OL residues, ONE phosphorus each) in two leaflets 38 A apart,
    # plus one cholesterol.
    for i in range(4):
        z = 19.0 if i < 2 else -19.0
        add("PC", [("P31", 30.97, [10.0 * i, 30.0, z])])
        add("PA", [("C12", 12.01, [10.0 * i, 31.0, z])])
        add("OL", [("C22", 12.01, [10.0 * i, 32.0, z])])
    add("CHL1", [("C1", 12.01, [5.0, 30.0, 19.0])])

    u = _universe(resnames, resids, names, masses, pos)
    u.dimensions = [60.0, 50.0, 80.0, 90.0, 90.0, 90.0]
    return u


def test_construct_residue_map_transfers_the_numbering(system, receptor_pdb):
    mapping = md.construct_residue_map(system, receptor_pdb)
    assert mapping.tolist() == list(range(69, 69 + len(RESNAMES)))


def test_construct_residue_map_refuses_a_mismatched_receptor(system, tmp_path):
    """A stale receptor.pdb must stop the analysis, not shift every residue number silently."""
    names, resids, masses, pos = _protein_atoms(["ALA", "ALA", "ALA"], 69)
    other = tmp_path / "wrong.pdb"
    _universe(["ALA", "ALA", "ALA"], resids, names, masses, pos, segid="R").atoms.write(str(other))
    with pytest.raises(ValueError, match="not the same receptor"):
        md.construct_residue_map(system, other)

    names, resids, masses, pos = _protein_atoms(["ALA"] * len(RESNAMES), 69)
    swapped = tmp_path / "swapped.pdb"
    _universe(["ALA"] * len(RESNAMES), resids, names, masses, pos,
              segid="R").atoms.write(str(swapped))
    with pytest.raises(ValueError, match="numbering transfer"):
        md.construct_residue_map(system, swapped)


def test_residue_min_distance_matches_a_direct_calculation(system):
    prot = system.select_atoms("protein")
    lig = system.select_atoms("resname LIG")
    dist = md.ResidueMinDistance(prot)(lig.positions, None)
    assert dist.size == len(RESNAMES)
    assert dist[2] == pytest.approx(3.0, abs=0.01)     # residue 3 sits 3.0 A from the ligand N
    assert dist.min() == pytest.approx(3.0, abs=0.01)
    # Same numbers the slow way.
    for i, res in enumerate(prot.residues):
        direct = np.linalg.norm(res.atoms.positions[:, None] - lig.positions[None], axis=2).min()
        assert dist[i] == pytest.approx(direct, abs=1e-4)


def test_occupancy_counts_frames_within_the_cutoff():
    d = np.array([[3.0, 9.0], [4.0, 9.0], [5.0, 9.0], [6.0, 9.0]])
    assert md.occupancy(d, cutoff=4.5).tolist() == [0.5, 0.0]


def test_water_bridges_finds_the_bridging_water(system):
    water_o = md.select_by_mass(system, *md.MASS_O, extra="resname WAT")
    lig_polar = system.select_atoms("resname LIG")
    prot = system.select_atoms("protein").residues
    polar = [r.atoms.positions for r in prot]
    counts = md.water_bridges(water_o.positions, polar, lig_polar.positions, None)
    assert counts[2] == 1        # construct 71 is bridged
    assert counts.sum() == 1     # and nothing else is


def test_lipid_counting_is_per_molecule_not_per_residue(system):
    n_phos, n_sterol = md.count_lipids(system)
    assert (n_phos, n_sterol) == (4, 1)   # 4 POPC = 12 residues, 1 cholesterol
    assert md.area_per_lipid((60.0, 50.0), n_phos + n_sterol) == pytest.approx(60 * 50 / 2.5)


def test_membrane_geometry(system):
    p = md.select_by_mass(system, *md.MASS_P).positions
    assert md.membrane_midplane(p) == pytest.approx(0.0, abs=1e-5)
    assert md.bilayer_thickness(p) == pytest.approx(38.0, abs=0.01)


def test_principal_components_of_a_one_dimensional_motion():
    """A single collective coordinate must land entirely in PC1."""
    n = 200
    coords = np.zeros((n, 4, 3))
    drive = np.linspace(-1.0, 1.0, n)
    coords[:, 0, 0] = drive
    coords[:, 1, 0] = -drive
    proj, var = md.principal_components(coords, n_components=2)
    assert proj.shape == (n, 2)
    assert var[0] == pytest.approx(1.0, abs=1e-6)
    assert abs(np.corrcoef(proj[:, 0], drive)[0, 1]) == pytest.approx(1.0, abs=1e-6)


def test_discover_builds_keeps_only_the_newest_per_system(tmp_path):
    for name in ("apo_ASP_20260901_120000", "apo_ASP_20260902_120000",
                 "ZH853_ASH_20260901_120000", "not_a_build", "ZH853_ASH_incomplete"):
        (tmp_path / name).mkdir()
        if "not_a_build" not in name and "incomplete" not in name:
            (tmp_path / name / "system.prmtop").touch()
    builds = md.discover_builds(tmp_path)
    assert [b.name for b in builds] == ["apo_ASP", "ZH853_ASH"]      # apo first, then ligands
    assert builds[0].stamp == "20260902_120000"                      # the newer of the two
    assert len(md.discover_builds(tmp_path, newest_only=False)) == 3


def test_build_dir_reads_its_system_json(tmp_path):
    path = tmp_path / "ZH853_ASP_20260901_120000"
    path.mkdir()
    (path / "system.prmtop").touch()
    (path / "system.json").write_text('{"ligand": "ZH853", "d250": "ASP", '
                                      '"ligand_resname": "LIG"}')
    (path / "prod_r2.dcd").touch()
    (path / "prod_r1.dcd").touch()
    build = md.discover_builds(tmp_path)[0]
    assert build.ligand_resname() == "LIG"
    assert [r for r, _ in build.replicas()] == ["prod_r1", "prod_r2"]

    apo = tmp_path / "apo_ASP_20260901_120000"
    apo.mkdir()
    (apo / "system.prmtop").touch()
    (apo / "system.json").write_text('{"ligand": "apo", "ligand_resname": ""}')
    assert md.discover_builds(tmp_path)[0].ligand_resname() is None


def test_read_state_log_parses_by_column_name(tmp_path):
    """Columns are found by name: a reordered or extended reporter must not shift the values."""
    log = tmp_path / "prod_r1.log"
    log.write_text(
        '#"Step","Time (ps)","Potential Energy (kJ/mole)","Temperature (K)",'
        '"Box Volume (nm^3)","Density (g/mL)","Speed (ns/day)"\n'
        "25000,100.0,-1.0e6,310.2,900.0,1.02,580.0\n"
        "50000,200.0,-1.1e6,309.8,901.0,1.03,585.0\n"
    )
    series = md.read_state_log(log)
    assert series["temperature"].tolist() == [310.2, 309.8]
    assert series["density"].tolist() == [1.02, 1.03]
    assert series["volume"].tolist() == [900.0, 901.0]
    assert series["speed"].tolist() == [580.0, 585.0]


def test_read_state_log_rejects_an_empty_log(tmp_path):
    empty = tmp_path / "empty.log"
    empty.write_text('#"Step","Temperature (K)"\n')
    with pytest.raises(ValueError, match="no data rows"):
        md.read_state_log(empty)


def test_protein_cross_section_and_corrected_area_per_lipid():
    """A 20x20 A square slab in the membrane slab: area 400, and it comes off the lipid area."""
    square = np.array([[0.0, 0.0, 0.0], [20.0, 0.0, 0.0], [20.0, 20.0, 0.0], [0.0, 20.0, 0.0],
                       [10.0, 10.0, 40.0]])  # the last point is outside the slab and ignored
    area = md.protein_cross_section(square, midplane_z=0.0)
    assert area == pytest.approx(400.0, abs=1e-6)
    gross = md.area_per_lipid((100.0, 100.0), 100)
    net = md.area_per_lipid((100.0, 100.0), 100, protein_area=400.0)
    assert gross == pytest.approx(200.0)
    assert net == pytest.approx(192.0)


def _write_reduced(root, system, replica, **series):
    """A minimal reduced-replica pair, as 01_reduce_trajectory.py writes them."""
    import json

    d = root / system
    d.mkdir(parents=True, exist_ok=True)
    ligand, d250 = system.rsplit("_", 1)
    (d / f"{replica}.json").write_text(json.dumps({
        "build": system, "ligand": ligand, "d250": d250, "replica": replica,
        "equilibration": {"t0_frames": 2}, "means": {"rmsd_ca_tm": {"mean": 1.5}}}))
    np.savez(d / f"{replica}.npz", **series)


def test_load_replicas_and_grouping(tmp_path):
    _write_reduced(tmp_path, "ZH853_ASH", "prod_r1", rmsd_ca_tm=np.arange(10.0))
    _write_reduced(tmp_path, "ZH853_ASH", "prod_r2", rmsd_ca_tm=np.arange(10.0))
    _write_reduced(tmp_path, "apo_ASH", "prod_r1", rmsd_ca_tm=np.arange(10.0))
    (tmp_path / "apo_ASH" / "orphan.json").write_text("{}")   # no .npz: ignored

    results = md.load_replicas(tmp_path)
    assert len(results) == 3
    systems = md.group_by_system(results)
    assert list(systems) == ["apo_ASH", "ZH853_ASH"]          # apo first, then ligands
    assert [r.replica for r in systems["ZH853_ASH"]] == ["prod_r1", "prod_r2"]

    rep = systems["ZH853_ASH"][0]
    assert rep.ligand == "ZH853" and rep.d250 == "ASH" and rep.t0 == 2
    assert rep.mean_of("rmsd_ca_tm") == 1.5
    assert np.isnan(rep.mean_of("not_measured"))
    assert rep.series("rmsd_ca_tm").tolist() == list(range(2, 10))        # t0 dropped
    assert rep.series("rmsd_ca_tm", equilibrated=False).size == 10
    with pytest.raises(KeyError):
        rep.series("nonexistent")

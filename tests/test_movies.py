"""Tests for the 02.12.00 movie exporter.

The exporter's job is to turn a wrapped, drifting, 5 GB DCD into a few MB that a viewer shows
CORRECTLY, and every way it can fail is silent: a torn lipid still renders, a ligand imaged a box
length away still renders, an atom written into two groups still renders. So the checks here are
on the things a picture cannot betray -- disjoint groups, whole molecules, drift removed by the
fit and left in without it -- on a synthetic system whose right answers are known by construction.

The step module is imported by path: `01_export_movie_trajectory.py` is a workflow script, not a
package module, and its name is not an identifier.
"""

from __future__ import annotations

import importlib.util
import sys

import MDAnalysis as mda
import numpy as np
import pytest

from zh853mor import md, paths

STEP = paths.SRC / "02.12.00_render_trajectory_movies" / "01_export_movie_trajectory.py"


def _load_step():
    spec = importlib.util.spec_from_file_location("export_movie_trajectory", STEP)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module          # dataclasses resolve types via sys.modules
    spec.loader.exec_module(module)
    return module


export = _load_step()

BOX = np.array([48.0, 48.0, 64.0])


def _system():
    """Protein (2 chains) + ligand + 3-residue lipids + OPC water + Na+, with real bonds.

    Deliberately mirrors the real build's awkward parts: a Lipid21 POPC spread over three bonded
    residues, an OPC water carrying a MASSLESS extra point, and monatomic ions.
    """
    names, types, masses, resnames, resids, atom_resindex, pos, bonds = [], [], [], [], [], [], [], []
    M = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "P": 30.974, "NA": 22.990, "EP": 0.0}

    def residue(resname, resid, atoms, xyz):
        ri = len(resnames)
        resnames.append(resname)
        resids.append(resid)
        first = len(names)
        for (name, elem), p in zip(atoms, xyz, strict=True):
            names.append(name)
            types.append(elem)
            masses.append(M[elem])
            atom_resindex.append(ri)
            pos.append(p)
        return first

    bb = [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")]
    off = np.array([[-1.2, 0, 0], [0, 0, 0], [1.2, 0, 0], [1.6, 0.9, 0]])
    # Chain 1: 24 residues climbing through the bilayer, construct-numbered from 69.
    first = None
    for k in range(24):
        f = residue("ALA", 69 + k, bb, np.array([2.0, 0.0, -20.0 + 1.7 * k]) + off)
        first = f if first is None else first
    bonds += [(i, i + 1) for i in range(first, len(names) - 1)]
    # Chain 2: a separate polypeptide below the membrane -- it must NOT be bonded to chain 1.
    first2 = None
    for k in range(6):
        f = residue("GLY", 1 + k, bb, np.array([-8.0, 0.0, -30.0 + 1.7 * k]) + off)
        first2 = f if first2 is None else first2
    bonds += [(i, i + 1) for i in range(first2, len(names) - 1)]

    ring = np.stack([2.0 * np.cos(np.linspace(0, 2 * np.pi, 8, endpoint=False)),
                     2.0 * np.sin(np.linspace(0, 2 * np.pi, 8, endpoint=False)),
                     np.zeros(8)], axis=1) + np.array([2.0, 0.0, 14.0])
    f = residue("LIG", 1, [(f"C{i}", "C") for i in range(8)], ring)
    bonds += [(f + i, f + i + 1) for i in range(7)] + [(f, f + 7)]

    lip = 1
    for leaflet, zhead in ((1, 15.0), (-1, -15.0)):
        for x in (-16.0, 0.0, 16.0):
            for y in (-16.0, 0.0, 16.0):
                head = np.array([x, y, zhead])
                f0 = residue("PC", lip, [("P", "P"), ("O1", "O"), ("N4", "N")],
                             head + np.array([[0, 0, 0], [1, 1, 0], [0, 0, 1.5 * leaflet]]))
                f1 = residue("OL", lip + 1, [(f"C{i}", "C") for i in range(4)],
                             head + np.array([[0.8, 0, -(2 + 2 * i) * leaflet] for i in range(4)]))
                # Named C4.. rather than D0..: read back from a PDB, MDAnalysis guesses the
                # element from the atom NAME, and "D0" guesses as deuterium (2.014) -- just heavy
                # enough to be ambiguous against the mass > 2.0 heavy-atom convention.
                f2 = residue("PA", lip + 2, [(f"C{i + 4}", "C") for i in range(4)],
                             head + np.array([[-0.8, 0, -(2 + 2 * i) * leaflet] for i in range(4)]))
                bonds += [(f0, f0 + 1), (f0, f0 + 2), (f0, f1), (f0, f2)]
                bonds += [(f1 + i, f1 + i + 1) for i in range(3)]
                bonds += [(f2 + i, f2 + i + 1) for i in range(3)]
                lip += 3

    rng = np.random.default_rng(3)
    waters = np.vstack([rng.uniform(-1, 1, (40, 3)) * np.array([20, 20, 28]),
                        np.array([[2.5, 0.5, 13.0], [1.0, 1.5, 15.5]])])   # two in the pocket
    for i, p in enumerate(waters):
        f = residue("WAT", i + 1, [("O", "O"), ("H1", "H"), ("H2", "H"), ("EPW", "EP")],
                    p + np.array([[0, 0, 0], [0.8, 0.6, 0], [-0.8, 0.6, 0], [0, 0.2, 0]]))
        bonds += [(f, f + 1), (f, f + 2), (f, f + 3)]
    for i, p in enumerate(rng.uniform(-1, 1, (6, 3)) * np.array([18, 18, 24])):
        residue("Na+", i + 1, [("Na+", "NA")], [p])

    n = len(names)
    u = mda.Universe.empty(n, n_residues=len(resnames), n_segments=1,
                           atom_resindex=np.array(atom_resindex),
                           residue_segindex=np.zeros(len(resnames), dtype=int), trajectory=True)
    u.add_TopologyAttr("name", names)
    u.add_TopologyAttr("type", types)
    u.add_TopologyAttr("mass", np.array(masses))
    u.add_TopologyAttr("resname", resnames)
    u.add_TopologyAttr("resid", np.array(resids))
    u.add_TopologyAttr("segid", ["SYSTEM"])
    u.add_bonds(bonds)
    u.atoms.positions = np.array(pos, dtype=np.float32)
    u.dimensions = [*BOX, 90.0, 90.0, 90.0]
    return u


def _wrap_by_molecule(positions, fragments, box):
    """Wrap each MOLECULE into [0, L) by its centroid -- what OpenMM writes into a DCD."""
    out = positions.copy()
    for frag in fragments:
        out[frag] -= np.floor(out[frag].mean(axis=0) / box) * box
    return out


DRIFT = np.array([9.0, 4.0, -6.0])     # planted whole-box translation, start to finish
N_FRAMES = 12


@pytest.fixture
def build(tmp_path):
    """A build directory indistinguishable from 02.10.00's, with a wrapped, drifting replica."""
    u = _system()
    path = tmp_path / "TST_ASP_20260101_000000"
    path.mkdir()
    base = u.atoms.positions.copy() + BOX / 2
    frags = [f.ix for f in u.atoms.fragments]
    rng = np.random.default_rng(11)

    with mda.Writer(str(path / "prod_r1.dcd"), n_atoms=len(u.atoms), dt=100.0) as W:
        for i in range(N_FRAMES):
            t = i / (N_FRAMES - 1)
            p = base + rng.normal(0, 0.05, base.shape) + DRIFT * t
            u.atoms.positions = _wrap_by_molecule(p, frags, BOX).astype(np.float32)
            u.trajectory.ts.dimensions = [*BOX, 90.0, 90.0, 90.0]
            W.write(u.atoms)

    # receptor.pdb is written in the OPM FRAME -- origin-centred, membrane midplane at z = 0 --
    # while the trajectory lives in the packed box running 0..L. Reproducing that offset is the
    # whole point of the fixture: a reference written in the box frame would make the pose RMSD
    # come out right even when the code compares the two frames against each other.
    u.atoms.positions = (base - BOX / 2).astype(np.float32)
    prot = u.select_atoms("protein")
    # receptor.pdb carries the protein AND the grafted ligand, as ligands.py --graft writes it:
    # it is the reference the ligand POSE RMSD is measured against, so a protein-only stand-in
    # would exercise only the fallback path.
    staged = prot + u.select_atoms("resname LIG")
    with mda.Writer(str(path / "receptor.pdb"), n_atoms=len(staged)) as W:
        W.write(staged)
    # MDAnalysis cannot write a prmtop; the file only has to EXIST for discover_builds, and the
    # readable topology is handed to the exporter with --topology, as on the smoke-test path.
    (path / "system.prmtop").touch()
    (path / "system.json").write_text('{"ligand_resname": "LIG"}\n')
    with mda.Writer(str(path / "system_top.pdb"), n_atoms=len(u.atoms), bonds="all") as W:
        W.write(u.atoms)
    return path


def _export(build, out, **kw):
    b = md.discover_builds(build.parent)[0]
    job = export.Job(b, "prod_r1", build / "prod_r1.dcd")
    return export.export_replica(job, out, topology=build / "system_top.pdb",
                                 frames=N_FRAMES, scan_frames=N_FRAMES, **kw)


def _movie(out):
    return mda.Universe(str(out / "prod_r1_movie.pdb"), str(out / "prod_r1_movie.xtc"))


def test_groups_are_disjoint_and_contiguous(build, tmp_path):
    """Every kept atom belongs to exactly one group, and each group is one slice of the file.

    Both renderers colour by those [start, stop) slices. An atom counted twice shifts every slice
    after it, and the movie is then coloured by an offset -- lipids drawn as water, and so on.
    """
    summary = _export(build, tmp_path / "out")
    groups = summary["atoms"]["groups"]
    ordered = [groups[g] for g in export.GROUP_ORDER]
    assert ordered[0][0] == 0
    for (_, stop), (start, _) in zip(ordered, ordered[1:], strict=False):
        assert stop == start                               # contiguous, no gaps, no overlap
    assert ordered[-1][1] == summary["atoms"]["total"]


def test_molecules_are_whole_in_every_frame(build, tmp_path):
    """A bond may stretch a little; it may never span the box, which is what a PBC tear looks like."""
    _export(build, tmp_path / "out")
    u = _movie(tmp_path / "out")
    pairs = np.array([[b[0].ix, b[1].ix] for b in u.bonds])
    assert len(pairs)                                       # CONECT records survived the write
    worst = 0.0
    for _ in u.trajectory:
        p = u.atoms.positions
        worst = max(worst, float(np.linalg.norm(p[pairs[:, 0]] - p[pairs[:, 1]], axis=1).max()))
    assert worst < 5.0                                      # << BOX.min() / 2 == 24 A


def test_superposition_removes_whole_box_drift(build, tmp_path):
    out = tmp_path / "out"
    _export(build, out)
    u = _movie(out)
    ca = u.select_atoms("protein and name CA")
    travel = []
    for _ in u.trajectory:
        travel.append(ca.positions.mean(axis=0))
    assert np.linalg.norm(np.array(travel)[-1] - np.array(travel)[0]) < 0.5


def test_no_align_keeps_the_drift_along_the_membrane_normal(build, tmp_path):
    """--no-align exists to show drift and tilt, so it must not centre along z as well.

    x and y are still centred: the box origin is arbitrary in the membrane plane, and a receptor
    that wanders out of frame is not a diagnosis.
    """
    out = tmp_path / "raw"
    summary = _export(build, out, align=False)
    assert summary["movie"]["align_selection"] == "none"
    u = _movie(out)
    ca = u.select_atoms("protein and name CA")
    z = []
    for _ in u.trajectory:
        z.append(ca.positions.mean(axis=0)[2])
    assert abs((z[-1] - z[0]) - DRIFT[2]) < 1.0


def test_lipids_come_in_whole_molecules(build, tmp_path):
    """Lipid21 splits a POPC over three residues; selecting by residue would keep a third of one."""
    _export(build, tmp_path / "out", lipid_cut=12.0)
    u = _movie(tmp_path / "out")
    lipids = u.select_atoms("chainID M")
    assert len(lipids)
    kept = {r.resname for r in lipids.residues}
    assert {"PC", "OL", "PA"} <= kept


def test_pocket_waters_are_ranked_not_taken_at_random(build, tmp_path):
    """The waters kept are the persistent ones, and their occupancy is reported, not implied."""
    summary = _export(build, tmp_path / "out", n_waters=2)
    waters = summary["waters"]
    assert len(waters) == 2
    assert all(w["occupancy"] > 0.5 for w in waters)
    assert waters == sorted(waters, key=lambda w: -w["occupancy"])


def test_topology_carries_chains_hetatm_and_construct_numbering(build, tmp_path):
    """Three things a bare coordinate dump gets wrong, each of which a viewer then shows wrongly."""
    summary = _export(build, tmp_path / "out")
    text = (tmp_path / "out" / "prod_r1_movie.pdb").read_text()
    assert summary["atoms"]["numbering"].startswith("human")

    chains = summary["atoms"]["chains"]
    assert chains["ligand"] == "L" and chains["water"] == "W"
    # The two polypeptides must be separate chains, or a cartoon is drawn from one into the other.
    protein_chains = set(chains["receptor"]) | set(chains["other_protein"])
    assert len(protein_chains) == 2
    assert not (protein_chains & {"L", "M", "W", "I"})

    lines = [ln for ln in text.splitlines() if ln.startswith(("ATOM", "HETATM"))]
    ligand = [ln for ln in lines if ln[21] == "L"]
    assert ligand and all(ln.startswith("HETATM") for ln in ligand)
    # Construct numbering, transferred from receptor.pdb: the receptor starts at 69, not tleap's 1.
    receptor = [int(ln[22:26]) for ln in lines
                if ln.startswith("ATOM") and ln[21] in set(chains["receptor"])]
    assert min(receptor) == 69


def test_the_receptor_chain_is_the_one_carrying_the_pocket(build, tmp_path):
    """Not the first chain in the file: tleap's order follows the packing input, not the biology."""
    summary = _export(build, tmp_path / "out")
    chains = summary["atoms"]["chains"]
    u = _movie(tmp_path / "out")
    receptor = u.select_atoms("chainID " + " ".join(chains["receptor"]))
    assert {int(r.resid) for r in receptor.residues} & set(md.KEY_RESIDUES)


def test_frame_spacing_is_reported(build, tmp_path):
    """The one number that says which events CANNOT be in the movie."""
    summary = _export(build, tmp_path / "out")
    movie = summary["movie"]
    assert movie["n_frames"] == N_FRAMES
    assert movie["frame_spacing_ns"] == pytest.approx(0.1, abs=1e-3)
    assert len(movie["time_ns"]) == N_FRAMES


# --- 02.11.00 ligand pose RMSD -------------------------------------------------------------------
# The reduction lives next door and is imported the same way. These cover the reference-frame bug
# that made lig_rmsd_pose report the ~90 A offset between the OPM frame receptor.pdb is written in
# and the 0..L simulation box, as a near-constant "pose RMSD" for the whole replica.

REDUCE = paths.SRC / "02.11.00_analyze_simulations" / "01_reduce_trajectory.py"


def _load_reduce():
    spec = importlib.util.spec_from_file_location("reduce_trajectory", REDUCE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _reduce(build, out, **kw):
    reduce_mod = _load_reduce()
    b = md.discover_builds(build.parent)[0]
    job = reduce_mod.Job(b, "prod_r1", build / "prod_r1.dcd")
    return reduce_mod, reduce_mod.reduce_replica(
        job, out, topology=build / "system_top.pdb", **kw)


def test_pose_rmsd_uses_the_deposited_reference_when_the_names_match(build, tmp_path):
    _, summary = _reduce(build, tmp_path / "out")
    assert summary["ligand_reference"] == "receptor.pdb (deposited pose)"
    pose = np.load(tmp_path / "out" / "prod_r1.npz")["lig_rmsd_pose"]
    assert np.isfinite(pose).all()
    assert pose[0] < 1.0                      # frame 0 IS the staged pose, up to the noise added


def test_pose_rmsd_falls_back_in_the_REFERENCE_frame_not_the_box_frame(build, tmp_path):
    """The regression. receptor.pdb is written in the OPM frame, centred near the origin; the
    trajectory is in a box running 0..L. A frame-0 fallback taken from the raw box coordinates is
    ~90 A from the aligned coordinate it gets compared against, and that offset -- not any motion
    of the ligand -- is what the series then reports, near-constant, for the whole replica.
    """
    # Break the name transfer the way the analog builds do: rename the reference ligand's atoms.
    ref = build / "receptor.pdb"
    lines = []
    for ln in ref.read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")) and ln[17:20].strip() == "LIG":
            ln = ln[:12] + f"{'X' + ln[12:16].strip():>4s}" + ln[16:]
        lines.append(ln)
    ref.write_text("\n".join(lines) + "\n")

    _, summary = _reduce(build, tmp_path / "out")
    assert summary["ligand_reference"] == "production frame 0 (receptor-aligned)"
    assert any("cannot transfer the deposited pose" in w for w in summary["warnings"])

    pose = np.load(tmp_path / "out" / "prod_r1.npz")["lig_rmsd_pose"]
    assert pose[0] == pytest.approx(0.0, abs=1e-6)    # frame 0 IS the fallback reference
    assert pose.max() < 25.0                          # motion, not a 90 A frame offset
    assert not any("FIRST frame" in w for w in summary["warnings"])


def test_a_reference_in_the_wrong_frame_is_caught_at_frame_zero(build, tmp_path):
    """The safety net: whatever the cause, a pose RMSD that STARTS tens of angstroms out is a
    broken reference, and saying so at frame 0 beats emitting a plausible-looking series."""
    reduce_mod = _load_reduce()
    ref = build / "receptor.pdb"
    # Move the reference ligand bodily, leaving the names intact so the transfer still succeeds.
    lines = []
    for ln in ref.read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")) and ln[17:20].strip() == "LIG":
            x = float(ln[30:38]) + 90.0
            ln = ln[:30] + f"{x:8.3f}" + ln[38:]
        lines.append(ln)
    ref.write_text("\n".join(lines) + "\n")

    b = md.discover_builds(build.parent)[0]
    job = reduce_mod.Job(b, "prod_r1", build / "prod_r1.dcd")
    summary = reduce_mod.reduce_replica(job, tmp_path / "out",
                                        topology=build / "system_top.pdb")
    assert summary["ligand_reference"] == "receptor.pdb (deposited pose)"
    assert any("FIRST frame" in w for w in summary["warnings"])


def test_duplicate_reference_atom_names_are_refused_not_silently_scrambled(build, tmp_path):
    """A repeated name kept whichever atom came last, giving a scrambled reference pose -- a wrong
    RMSD that still looks like an RMSD, which is worse than no RMSD."""
    ref = build / "receptor.pdb"
    lines, first = [], True
    for ln in ref.read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")) and ln[17:20].strip() == "LIG" and not first:
            ln = ln[:12] + f"{'C0':>4s}" + ln[16:]        # collide every atom onto one name
        elif ln.startswith(("ATOM", "HETATM")) and ln[17:20].strip() == "LIG":
            first = False
        lines.append(ln)
    ref.write_text("\n".join(lines) + "\n")

    _, summary = _reduce(build, tmp_path / "out")
    assert summary["ligand_reference"] == "production frame 0 (receptor-aligned)"

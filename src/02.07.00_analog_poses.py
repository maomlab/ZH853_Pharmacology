#!/usr/bin/env python3
"""3-D poses for the three ZH analogs by constrained embedding on the common scaffold.

ZH853 is the only one of the four with an experimental pose. ZH850/ZH831/ZH809 are close analogs
of it -- same Tyr-cyclo[...] macrocyclic scaffold, differing by a truncation (ZH850 drops the
Gly-NH2 tail), a residue swap (ZH831 Trp->Phe with the cycle reversed), or one CH2 in the ring
(ZH809 Glu->Asp). For that degree of similarity the defensible way to place them is to inherit the
crystallographic binding mode rather than re-dock: find the maximum common substructure with
ZH853, pin those atoms to the ZH853 coordinates, and embed + minimise only what is left.

Everything happens in the OPM-oriented frame, using the ligand out of
`02.05.00_oriented/complex_oriented.pdb` as the template, so the emitted poses are already in the
frame the membrane build expects. The SDF frame itself is irrelevant to parameterization
(antechamber needs topology and a conformer; a rigid transform changes neither the charges nor the
internal geometry), but keeping one frame throughout means the pose and the complex cannot drift
apart.

Outputs, which is exactly what src/02.10.00_slurm_bundle/ligands.py looks for:
    intermediate/02.04.00_ligand/<NAME>_prepared.sdf
    intermediate/02.05.00_oriented/complex_<NAME>_oriented.pdb

Run: ``python src/02.07.00_analog_poses.py``  (or ``make prep-analogs``).
"""

from __future__ import annotations

import sys
from datetime import date

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import AllChem, rdFMCS  # noqa: E402
from rdkit.Chem import rdMolAlign  # noqa: E402

from zh853mor import chem, paths  # noqa: E402

RDLogger.DisableLog("rdApp.*")

TEMPLATE = "ZH853"
ANALOGS = ["ZH850", "ZH831", "ZH809"]
LIG_RESNAME = "L01"          # what the deposited pose calls it; fix_ligand.py renames it to LIG
N_CONFS = 25                 # embedding attempts; the best-restrained one is kept
CLASH_A = 2.0                # heavy-atom contact closer than this to the receptor is a clash
MCS_TIMEOUT_S = 60
BOND_RANGE_A = (0.9, 1.8)    # plausible heavy-heavy bond length; see prune_mapping()
TETHER_K = 50.0              # kcal/mol/A^2 spring holding the scaffold during stage-2 relaxation


def protonate(mol: Chem.Mol) -> Chem.Mol:
    """+1 at the Tyr1 alpha-amine, the group that salt-bridges D3.32 (Asp149).

    Same rule as 02.04.00_ligand_prep.py: a primary aliphatic N that is not an amide.
    """
    amine = Chem.MolFromSmarts("[NX3;H2;!$(NC=O)]")
    for (idx,) in mol.GetSubstructMatches(amine):
        mol.GetAtomWithIdx(idx).SetFormalCharge(1)
    Chem.SanitizeMol(mol)
    return mol


def template_from_complex(path) -> Chem.Mol:
    """ZH853 heavy atoms out of the oriented complex, with bond orders from its SMILES."""
    lines = [l for l in path.read_text().splitlines()
             if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == LIG_RESNAME]
    if not lines:
        raise SystemExit(f"ERROR: no {LIG_RESNAME} residue in {path}")
    pdb_mol = Chem.MolFromPDBBlock("\n".join(lines) + "\nEND\n",
                                   sanitize=False, proximityBonding=True)
    ref = Chem.MolFromSmiles(chem.ANALOGS[TEMPLATE][0])
    mol = AllChem.AssignBondOrdersFromTemplate(ref, pdb_mol)
    return protonate(mol)


def mcs_map(analog: Chem.Mol, template: Chem.Mol):
    """(analog_idx, template_idx, n_mcs) for the maximum common substructure, or None.

    completeRingsOnly is off deliberately: ZH809's macrocycle is one CH2 smaller than ZH853's, so
    requiring whole rings to match would discard the entire scaffold -- the very atoms whose
    placement is known -- and leave nothing to constrain.
    """
    res = rdFMCS.FindMCS(
        [analog, template],
        ringMatchesRingOnly=True, completeRingsOnly=False,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
        timeout=MCS_TIMEOUT_S,
    )
    if res.canceled or res.numAtoms < 8:
        return None
    patt = Chem.MolFromSmarts(res.smartsString)
    a_match, t_match = analog.GetSubstructMatch(patt), template.GetSubstructMatch(patt)
    if not a_match or not t_match:
        return None
    return list(a_match), list(t_match), res.numAtoms


def prune_mapping(analog, template, a_match, t_match):
    """Drop scaffold atoms whose inherited coordinates would violate bond geometry.

    An MCS match is a TOPOLOGICAL statement; it does not promise that atoms bonded in the analog
    map to atoms that are adjacent in the template. Where the macrocycle is traversed differently
    -- ZH831 reverses the cycle -- a valid match can put a bonded pair 4.5 A apart. Inheriting
    those coordinates hands the force field an impossible geometry (MMFF came back at
    7.5e4 kcal/mol for ZH831 before this check existed). So keep only the subset whose induced
    bond lengths are physical, and rebuild the rest properly.
    """
    import numpy as np
    tconf = template.GetConformer()
    pos = {a: np.array([tconf.GetAtomPosition(t).x, tconf.GetAtomPosition(t).y,
                        tconf.GetAtomPosition(t).z]) for a, t in zip(a_match, t_match)}
    mapping = dict(zip(a_match, t_match))
    lo, hi = BOND_RANGE_A
    while True:
        bad = {}
        for b in analog.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if i in mapping and j in mapping:
                d = float(np.linalg.norm(pos[i] - pos[j]))
                if not (lo <= d <= hi):
                    bad[i] = bad.get(i, 0) + 1
                    bad[j] = bad.get(j, 0) + 1
        if not bad:
            break
        # Remove the atom implicated in the most bad bonds (ties -> highest index, deterministic).
        worst = max(sorted(bad), key=lambda k: bad[k])
        del mapping[worst]
    a_keep = sorted(mapping)
    return a_keep, [mapping[a] for a in a_keep]


def embed(analog, template, a_match, t_match):
    """Build the analog pose by inheriting the ZH853 scaffold coordinates outright.

    The MCS covers 49-55 of each analog's 52-55 heavy atoms -- ZH850 and ZH809 are entirely
    contained in it -- so distance geometry is the wrong instrument: constraining ~all atoms makes
    the embedding either fail or fight itself (ConstrainedEmbed failed outright on all three, and a
    coordMap embed then blew up MMFF's line search). With overlap this high the scaffold is simply
    COPIED, exactly, and only the handful of genuinely new atoms is rebuilt.

    Returns (mol with H, n_rebuilt, mmff_energy).
    """
    import numpy as np
    tconf = template.GetConformer()
    conf = Chem.Conformer(analog.GetNumAtoms())
    placed = {}
    for a, t in zip(a_match, t_match):
        pos = tconf.GetAtomPosition(t)
        conf.SetAtomPosition(a, pos)
        placed[a] = np.array([pos.x, pos.y, pos.z])

    # Seed each new atom just off a placed neighbour; the restrained minimisation below is what
    # actually positions it. Deterministic offsets so a rerun reproduces the pose.
    rng = np.random.default_rng(0)
    todo = [a.GetIdx() for a in analog.GetAtoms() if a.GetIdx() not in placed]
    for _ in range(len(todo) + 1):
        if not todo:
            break
        for idx in list(todo):
            nbrs = [n.GetIdx() for n in analog.GetAtomWithIdx(idx).GetNeighbors() if n.GetIdx() in placed]
            if not nbrs:
                continue
            v = rng.normal(size=3)
            p = placed[nbrs[0]] + 1.5 * v / np.linalg.norm(v)
            conf.SetAtomPosition(idx, Chem.rdGeometry.Point3D(*p))
            placed[idx] = p
            todo.remove(idx)
    if todo:
        raise RuntimeError(f"{len(todo)} atom(s) are not connected to the common scaffold")

    mol = Chem.Mol(analog)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    molh = Chem.AddHs(mol, addCoords=True)   # AddHs appends, so heavy indices are preserved

    # Relax hydrogens and the rebuilt atoms with the inherited scaffold held rigid, so the
    # crystallographic binding mode is preserved exactly rather than approximately.
    energy = float("nan")
    props = AllChem.MMFFGetMoleculeProperties(molh)
    ff = AllChem.MMFFGetMoleculeForceField(molh, props) if props is not None else \
        AllChem.UFFGetMoleculeForceField(molh)
    if ff is not None:
        for a in a_match:
            ff.AddFixedPoint(a)
        ff.Initialize()
        try:
            ff.Minimize(maxIts=2000)
            energy = ff.CalcEnergy()
        except RuntimeError as exc:      # BFGS line-search failures are not fatal here: the
            print(f"  note: minimisation stopped early ({exc.__class__.__name__}); "
                  "the inherited scaffold is unaffected")
    # Stage 2: let residual strain relax with the scaffold TETHERED rather than frozen. Stage 1
    # can leave a rebuilt atom fighting a bond it cannot reach while its neighbours are immovable
    # (ZH831, whose reversed cycle makes it the least ZH853-like of the three). Springs let the
    # geometry heal by a few tenths of an angstrom while the binding mode is still held.
    if ff is not None:
        ff2 = AllChem.MMFFGetMoleculeForceField(molh, props) if props is not None else \
            AllChem.UFFGetMoleculeForceField(molh)
        tconf_ = template.GetConformer()
        for a, t in zip(a_match, t_match):
            pos = tconf_.GetAtomPosition(t)
            anchor = ff2.AddExtraPoint(pos.x, pos.y, pos.z, fixed=True) - 1
            ff2.AddDistanceConstraint(anchor, a, 0, 0, TETHER_K)
        ff2.Initialize()
        try:
            ff2.Minimize(maxIts=2000)
            energy = ff2.CalcEnergy()
        except RuntimeError:
            pass
    drift = rdMolAlign.CalcRMS(molh, template, map=[list(zip(a_match, t_match))]) \
        if len(a_match) >= 3 else float("nan")
    return molh, len(analog.GetAtoms()) - len(a_match), energy, drift


def strain_baseline(template) -> float:
    """MMFF energy of the ZH853 template pose itself, H relaxed -- the only meaningful yardstick.

    Absolute MMFF energies are not interpretable across molecules, but the analogs differ from the
    template by a few atoms, so 'how much worse than the pose we trust' is a fair comparison.
    """
    molh = Chem.AddHs(template, addCoords=True)
    props = AllChem.MMFFGetMoleculeProperties(molh)
    if props is None:
        return float("nan")
    ff = AllChem.MMFFGetMoleculeForceField(molh, props)
    for a in range(template.GetNumAtoms()):
        ff.AddFixedPoint(a)
    ff.Initialize()
    try:
        ff.Minimize(maxIts=2000)
    except RuntimeError:
        pass
    return ff.CalcEnergy()


def receptor_heavy(path):
    import numpy as np
    return np.array([[float(l[30 + 8 * i:38 + 8 * i]) for i in range(3)]
                     for l in path.read_text().splitlines()
                     if l.startswith("ATOM") and (l[76:78].strip() or l[12:16].strip()[0]) != "H"])


def write_outputs(name, molh, oriented_complex, lig_dir, ori_dir) -> tuple[int, float]:
    """Write <name>_prepared.sdf and complex_<name>_oriented.pdb; return (n_heavy, min clash)."""
    import numpy as np
    sdf = lig_dir / f"{name}_prepared.sdf"
    molh.SetProp("_Name", name)
    w = Chem.SDWriter(str(sdf)); w.write(molh); w.close()

    conf = molh.GetConformer()
    heavy = [a.GetIdx() for a in molh.GetAtoms() if a.GetAtomicNum() != 1]
    # Receptor block copied verbatim from the ZH853 oriented complex: the bundle checks that the
    # complex and receptorR_oriented.pdb share a CA frame, and copying guarantees they do.
    rec = [l for l in oriented_complex.read_text().splitlines()
           if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() != LIG_RESNAME]
    lig_lines, coords = [], []
    for i, idx in enumerate(heavy, 1):
        p = conf.GetAtomPosition(idx)
        el = molh.GetAtomWithIdx(idx).GetSymbol()
        # Heavy-atom ORDER here must match the SDF, because fix_ligand.py maps positionally.
        lig_lines.append(f"HETATM{i:>5} {el + str(i):<4} {LIG_RESNAME} E   1    "
                         f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00          {el:>2}")
        coords.append([p.x, p.y, p.z])

    out = ori_dir / f"complex_{name}_oriented.pdb"
    out.write_text("\n".join(rec + lig_lines) + "\nEND\n")

    d = np.linalg.norm(np.array(coords)[:, None, :] - receptor_heavy(out)[None, :, :], axis=2)
    return len(heavy), float(d.min())


def main() -> int:
    lig_dir = paths.ensure_dir(paths.INTERMEDIATE / "02.04.00_ligand")
    ori_dir = paths.INTERMEDIATE / "02.05.00_oriented"
    cx_path = ori_dir / "complex_oriented.pdb"
    if not cx_path.exists():
        print(f"ERROR: {cx_path} not found.", file=sys.stderr)
        print("  This is the only input this step needs (that and rdkit). Two ways to get it:",
              file=sys.stderr)
        print("   - in the LOCAL analysis env:  make prep-receptor prep-orient", file=sys.stderr)
        print("     (needs openmm + pdbfixer, which the cluster prep env does not carry by design)",
              file=sys.stderr)
        print("   - or copy the file across; intermediate/ is gitignored, so a git pull will not",
              file=sys.stderr)
        print("     bring it:  scp <local>/intermediate/02.05.00_oriented/complex_oriented.pdb \\",
              file=sys.stderr)
        print(f"                    <host>:{cx_path}", file=sys.stderr)
        return 1

    template = template_from_complex(cx_path)
    baseline = strain_baseline(template)
    print(f"template {TEMPLATE}: {template.GetNumAtoms()} heavy atoms from {cx_path.name}, "
          f"MMFF {baseline:.1f} kcal/mol (the yardstick for the analogs below)\n")

    rows, failed = [], 0
    for name in ANALOGS:
        mol = protonate(Chem.MolFromSmiles(chem.ANALOGS[name][0]))
        got = mcs_map(mol, template)
        if got is None:
            print(f"{name}: FAILED -- no usable MCS with {TEMPLATE}", file=sys.stderr)
            failed += 1
            continue
        a_match, t_match, n_mcs = got
        a_match, t_match = prune_mapping(mol, template, a_match, t_match)
        n_core = len(a_match)
        if n_core < n_mcs:
            print(f"{name}: pruned {n_mcs - n_core} MCS atom(s) whose inherited coordinates "
                  "violated bond geometry (the cycle is traversed differently); rebuilding them")
        if n_core < 8:
            print(f"{name}: FAILED -- only {n_core} atoms survive as a usable scaffold; "
                  "this analog needs docking, not scaffold transfer", file=sys.stderr)
            failed += 1
            continue
        try:
            molh, n_new, energy, drift = embed(mol, template, a_match, t_match)
        except RuntimeError as exc:
            print(f"{name}: FAILED -- {exc}", file=sys.stderr)
            failed += 1
            continue
        n_heavy, clash = write_outputs(name, molh, cx_path, lig_dir, ori_dir)
        frac = 100.0 * n_core / mol.GetNumAtoms()
        flag = "" if clash >= CLASH_A else f"  <-- CLASH (< {CLASH_A} A)"
        print(f"{name}: scaffold {n_core}/{mol.GetNumAtoms()} ({frac:.0f}%) inherited, "
              f"{n_new} rebuilt | MMFF {energy:8.1f} ({energy - baseline:+.0f} vs template) | "
              f"scaffold drift {drift:.2f} A | closest receptor contact {clash:.2f} A{flag}")
        rows.append((name, n_core, mol.GetNumAtoms(), frac, n_new, n_heavy, clash))

    print()
    for name, n_core, n_tot, frac, n_new, n_heavy, clash in rows:
        print(f"  wrote {name}_prepared.sdf + complex_{name}_oriented.pdb")
    if failed:
        print(f"\n{failed} analog(s) failed.", file=sys.stderr)
        return 1

    bad = [r for r in rows if r[6] < CLASH_A]
    if bad:
        print(f"\nWARNING: {len(bad)} pose(s) clash with the receptor (< {CLASH_A} A). "
              "Equilibration may not resolve a buried overlap -- inspect before building.")
    print(f"\nNext: parameterize and build, e.g.")
    print(f"    ./ligand_resp/run_resp.sh ZH850 && LIGAND=ZH850 ./01_build_system.sh")
    return 0


if __name__ == "__main__":
    sys.exit(main())

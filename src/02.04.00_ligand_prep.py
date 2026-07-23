#!/usr/bin/env python3
"""Prepare ZH853 for force-field parameterization (Phase 2).

Takes the deposited ligand coordinates, transfers bond orders/aromaticity from the reference
SMILES, protonates to the physiological +1 state (the Tyr1 alpha-amine that salt-bridges
D149), adds hydrogens on the 3D pose, and writes an SDF ready for antechamber/OpenFF. Also
emits the parameterization commands for both routes (GAFF2/AM1-BCC quick; RESP rigorous).

Run: ``python src/02.04.00_ligand_prep.py``  (or ``make prep-ligand``).
"""

from __future__ import annotations

import sys
from datetime import date

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402

from zh853mor import chem, comparators, paths  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def build_protonated_ligand() -> tuple[Chem.Mol, int]:
    """Return (ligand mol with 3D coords, explicit H, +1 charge; net_charge)."""
    cx = comparators.load_complex("ZH853")
    lig = cx.ligand
    # PDB block for the 59 heavy atoms (perceive connectivity from geometry)
    block_lines = []
    for i, a in enumerate(lig.atoms):
        x, y, z = a.position
        el = a.name[0]
        block_lines.append(
            f"HETATM{i + 1:>5} {a.name:<4} L01 E   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2}"
        )
    pdb_mol = Chem.MolFromPDBBlock("\n".join(block_lines) + "\nEND\n",
                                   sanitize=False, proximityBonding=True)
    template = Chem.MolFromSmiles(chem.ANALOGS["ZH853"][0])
    mol = AllChem.AssignBondOrdersFromTemplate(template, pdb_mol)

    # protonate the primary aliphatic amine (Tyr1 alpha-amine; not an amide N) -> +1
    amine = Chem.MolFromSmarts("[NX3;H2;!$(NC=O)]")
    matches = mol.GetSubstructMatches(amine)
    for (n_idx,) in matches:
        mol.GetAtomWithIdx(n_idx).SetFormalCharge(1)
    Chem.SanitizeMol(mol)
    molh = Chem.AddHs(mol, addCoords=True)
    return molh, Chem.GetFormalCharge(molh)


def main() -> int:
    today = f"{date.today():%Y%m%d}"
    work = paths.ensure_dir(paths.INTERMEDIATE / "02.04.00_ligand")
    molh, net_charge = build_protonated_ligand()

    sdf = work / "ZH853_prepared.sdf"
    writer = Chem.SDWriter(str(sdf))
    molh.SetProp("_Name", "ZH853")
    writer.write(molh)
    writer.close()
    pdb = work / "ZH853_prepared.pdb"
    Chem.MolToPDBFile(molh, str(pdb))

    n_h = sum(a.GetAtomicNum() == 1 for a in molh.GetAtoms())
    lines = [f"# ZH853 ligand preparation — {today}", "",
             f"- Bond orders/aromaticity transferred from the reference SMILES (template match OK).",
             f"- Protonation: **net charge {net_charge:+d}** (Tyr1 alpha-amine protonated; this is the "
             "group that salt-bridges D149). No free carboxylate — the Glu is in the macrocyclic lactam.",
             f"- Hydrogens added on the 3D pose: {n_h} H ({molh.GetNumAtoms()} atoms total).",
             "", "## Outputs (intermediate/02.04.00_ligand/)",
             f"- `{sdf.name}` — parameterization input (SDF, explicit H, correct bond orders, +1).",
             f"- `{pdb.name}` — same, PDB.",
             "", "## Parameterization routes (run on the cluster; see the SLURM bundle)",
             "**A. Quick (GAFF2 + AM1-BCC)** — good for equilibration/MD, minutes:",
             "```bash",
             "antechamber -i ZH853_prepared.sdf -fi sdf -o ZH853.mol2 -fo mol2 \\",
             "            -c bcc -nc 1 -at gaff2 -rn LIG",
             "parmchk2 -i ZH853.mol2 -f mol2 -o ZH853.frcmod -s gaff2",
             "```",
             "**B. Rigorous (multi-conformer RESP, HF/6-31G(d))** — for FEP, matches Amber charge "
             "derivation (SPECIFICATION D-3):",
             "```bash",
             "# generate conformers -> Psi4/Gaussian ESP -> RESP fit (2-stage) -> mol2 with RESP charges",
             "# (script staged in the SLURM bundle: src/02.10.00_slurm_bundle/ligand_resp/)",
             "```",
             "**C. OpenFF cross-check** — independent SMIRNOFF parameters via openmmforcefields "
             "SystemGenerator (`SMIRNOFFTemplateGenerator`, openff-2.x).",
             "", "> The macrocyclic backbone + D-amino acids are handled as one small molecule here "
             "(59 atoms, +1). For the residue-library alternative (capped-fragment RESP per non-canonical "
             "residue), see the bundle README. Validate cis/trans amide and ring conformers with enhanced "
             "sampling (D-5 risk).",
             ]
    (paths.PRODUCT / f"02.04.00_ligand_prep_{today}.md").write_text("\n".join(lines) + "\n")
    print(f"Prepared ZH853: net charge {net_charge:+d}, {n_h} H added.")
    print(f"Outputs in {work}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

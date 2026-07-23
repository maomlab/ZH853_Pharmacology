"""Cheminformatics for ZH853 and analogs (Objective 3: drug-likeness).

SMILES are taken from OBJECTIVES.md and treated as authoritative (the analog *names* there
have Trp/Phe inconsistencies — see SPECIFICATION OQ-5/D-6). Provides 2D and best-effort 3D
descriptors, beyond-rule-of-5 classification, and builders for the modification strategies
(N-methylation, lipidation, halogenation).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, Descriptors3D, rdFreeSASA, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

# name -> (SMILES, description). ZH853 is the deposited ligand; the other three are the
# Objective-4 FEP analogs. All are Tyr-cyclo[...]-NH2 endomorphin-1 macrocycles.
ANALOGS: dict[str, tuple[str, str]] = {
    "ZH853": (
        "NC(=O)CNC(=O)[C@@H]5CCC(=O)NCCCC[C@@H](NC(=O)[C@@H](N)Cc1ccc(O)cc1)C(=O)"
        "N[C@@H](Cc2c[nH]c3ccccc23)C(=O)N[C@@H](Cc4ccccc4)C(=O)N5",
        "Tyr-cyclo[D-Lys-Trp-Phe-Glu]-Gly-NH2 (deposited ligand)",
    ),
    "ZH850": (
        "NC([C@@H]1CCC(NCCCC[C@H](C(N[C@H](C(N[C@H](C(N1)=O)Cc2ccccc2)=O)"
        "Cc(c[nH]3)c4c3cccc4)=O)NC([C@H](Cc5ccc(O)cc5)N)=O)=O)=O",
        "Tyr-cyclo[D-Lys-Trp-Phe-Glu]-NH2 (analog 1)",
    ),
    "ZH831": (
        "NC(=O)[C@@H]4CCCCNC(=O)CC[C@@H](NC(=O)[C@@H](N)Cc1ccc(O)cc1)C(=O)"
        "N[C@@H](Cc2ccccc2)C(=O)N[C@@H](Cc3ccccc3)C(=O)N4",
        "Tyr-cyclo[D-Glu-Phe-Phe-Lys]-NH2 (analog 2)",
    ),
    "ZH809": (
        "NC(=O)[C@@H]5CC(=O)NCCCC[C@@H](NC(=O)[C@@H](N)Cc1ccc(O)cc1)C(=O)"
        "N[C@@H](Cc2c[nH]c3ccccc23)C(=O)N[C@@H](Cc4ccccc4)C(=O)N5",
        "Tyr-cyclo[D-Lys-Trp-Phe-Asp]-NH2 (analog 3)",
    ),
}

# Lipidation linkers (dummy atom [*] marks the attachment bond to a ligand amine).
LINKERS: dict[str, tuple[str, str]] = {
    "palmitoyl": ("[*]C(=O)CCCCCCCCCCCCCCC", "C16 fatty acyl (simple lipidation)"),
    "semaglutide-like": (
        "[*]C(=O)CC[C@@H](C(=O)O)NC(=O)COCCOCCNC(=O)COCCOCCNC(=O)CCCCCCCCCCCCCCCCC(=O)O",
        "gamma-Glu-2xAEEA-C18-diacid (albumin-binding, GLP-1 style half-life extension)",
    ),
}


@dataclass
class Props:
    """2D and (best-effort) 3D physicochemical descriptors."""

    name: str
    smiles: str
    mw: float = 0.0
    tpsa: float = 0.0
    clogp: float = 0.0
    hbd: int = 0
    hba: int = 0
    rotb: int = 0
    aromatic_rings: int = 0
    fsp3: float = 0.0
    formal_charge: int = 0
    # 3D (may be NaN if embedding fails)
    rgyr: float = float("nan")
    psa3d_frac: float = float("nan")
    intramol_hbonds: float = float("nan")
    bro5_flags: list[str] = field(default_factory=list)


def mol(name_or_smiles: str) -> Chem.Mol:
    """Parse a SMILES (or a known analog name) to an RDKit mol."""
    smi = ANALOGS[name_or_smiles][0] if name_or_smiles in ANALOGS else name_or_smiles
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise ValueError(f"Could not parse SMILES: {name_or_smiles}")
    return m


def descriptors_2d(name: str, m: Chem.Mol) -> Props:
    p = Props(name=name, smiles=Chem.MolToSmiles(m))
    p.mw = Descriptors.MolWt(m)
    p.tpsa = rdMolDescriptors.CalcTPSA(m)
    p.clogp = Descriptors.MolLogP(m)
    p.hbd = rdMolDescriptors.CalcNumHBD(m)
    p.hba = rdMolDescriptors.CalcNumHBA(m)
    p.rotb = rdMolDescriptors.CalcNumRotatableBonds(m)
    p.aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(m)
    p.fsp3 = rdMolDescriptors.CalcFractionCSP3(m)
    p.formal_charge = Chem.GetFormalCharge(m)
    p.bro5_flags = _bro5(p)
    return p


def add_descriptors_3d(p: Props, m: Chem.Mol, n_conf: int = 5, seed: int = 0xF00D) -> Props:
    """Embed a few macrocycle-aware conformers; report the lowest-energy one's 3D descriptors."""
    mh = Chem.AddHs(m)
    params = AllChem.ETKDGv3()
    params.useMacrocycleTorsions = True
    params.randomSeed = seed
    cids = AllChem.EmbedMultipleConfs(mh, numConfs=n_conf, params=params)
    if not len(cids):
        return p
    energies = []
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mh, AllChem.MMFFGetMoleculeProperties(mh), confId=cid)
        if ff is None:
            continue
        ff.Minimize()
        energies.append((ff.CalcEnergy(), cid))
    if not energies:
        return p
    best = min(energies)[1]
    p.rgyr = float(Descriptors3D.RadiusOfGyration(mh, confId=best))
    p.psa3d_frac = _polar_sasa_fraction(mh, best)
    p.intramol_hbonds = _intramolecular_hbonds(mh, best)
    return p


def _polar_sasa_fraction(mh: Chem.Mol, conf_id: int) -> float:
    radii = rdFreeSASA.classifyAtoms(mh)
    total = rdFreeSASA.CalcSASA(mh, radii, confIdx=conf_id)
    if total <= 0:
        return float("nan")
    polar = 0.0
    for a in mh.GetAtoms():
        is_polar = a.GetSymbol() in ("N", "O") or (
            a.GetSymbol() == "H" and any(nb.GetSymbol() in ("N", "O") for nb in a.GetNeighbors())
        )
        if is_polar:
            polar += a.GetPropsAsDict().get("SASA", 0.0)
    return polar / total


def _intramolecular_hbonds(mh: Chem.Mol, conf_id: int, cutoff: float = 2.5) -> int:
    conf = mh.GetConformer(conf_id)
    pos = conf.GetPositions()
    donor_h = [
        (a.GetIdx(), a.GetNeighbors()[0].GetIdx())
        for a in mh.GetAtoms()
        if a.GetSymbol() == "H" and a.GetNeighbors() and a.GetNeighbors()[0].GetSymbol() in ("N", "O")
    ]
    acceptors = [a.GetIdx() for a in mh.GetAtoms() if a.GetSymbol() in ("N", "O")]
    n = 0
    for h, donor in donor_h:
        for acc in acceptors:
            if acc != donor and np.linalg.norm(pos[h] - pos[acc]) < cutoff:
                n += 1
                break
    return n


def _bro5(p: Props) -> list[str]:
    """Flag Lipinski / beyond-rule-of-5 violations (Doak bRo5 thresholds noted in report)."""
    flags = []
    if p.mw > 500:
        flags.append("MW>500")
    if p.tpsa > 140:
        flags.append("TPSA>140")
    if p.hbd > 5:
        flags.append("HBD>5")
    if p.hba > 10:
        flags.append("HBA>10")
    if p.clogp > 5:
        flags.append("cLogP>5")
    return flags


# ---- modification builders ---------------------------------------------------

def n_methyl_sites(m: Chem.Mol) -> list[int]:
    """Atom indices of backbone secondary-amide nitrogens available for N-methylation."""
    patt = Chem.MolFromSmarts("[NX3;H1][CX3]=[OX1]")
    return [match[0] for match in m.GetSubstructMatches(patt)]


def add_n_methyls(m: Chem.Mol, n_indices: list[int]) -> Chem.Mol:
    """Return a copy with a methyl added to each given amide nitrogen."""
    rw = Chem.RWMol(m)
    for idx in n_indices:
        c = rw.AddAtom(Chem.Atom(6))
        rw.AddBond(idx, c, Chem.BondType.SINGLE)
        rw.GetAtomWithIdx(idx).SetNumExplicitHs(0)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out


def conjugate_at_primary_amide(m: Chem.Mol, linker_smiles: str) -> Chem.Mol:
    """Attach a linker (with a [*] dummy) to a primary amide nitrogen of ``m``."""
    patt = Chem.MolFromSmarts("[NX3;H2][CX3]=[OX1]")
    matches = m.GetSubstructMatches(patt)
    if not matches:
        raise ValueError("no primary amide site for conjugation")
    n_idx = matches[0][0]
    frag = Chem.MolFromSmiles(linker_smiles)
    dummy = next(a.GetIdx() for a in frag.GetAtoms() if a.GetAtomicNum() == 0)
    attach = frag.GetAtomWithIdx(dummy).GetNeighbors()[0].GetIdx()
    frag = Chem.DeleteSubstructs(frag, Chem.MolFromSmarts("[#0]"))
    attach -= 1 if attach > dummy else 0  # index shift after dummy removal

    combo = Chem.RWMol(Chem.CombineMols(m, frag))
    offset = m.GetNumAtoms()
    combo.AddBond(n_idx, offset + attach, Chem.BondType.SINGLE)
    combo.GetAtomWithIdx(n_idx).SetNumExplicitHs(1)
    out = combo.GetMol()
    Chem.SanitizeMol(out)
    return out

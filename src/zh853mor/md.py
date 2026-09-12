"""Trajectory-analysis primitives for the production membrane simulations.

Shared by `src/02.11.00_analyze_simulations/`: the per-replica reduction runs on the cluster
(where the trajectories are), the aggregation and figures run locally on what it writes. Both
import this module, so a selection or a cutoff is defined once.

Three things here are easy to get wrong and are therefore done centrally:

**Residue numbering.** tleap renumbers the system from 1, so a prmtop residue id is NOT the human
OPRM1 number every other file in this project speaks (D149, E231, H299 ...). The staged
``receptor.pdb`` in each build directory keeps the construct numbering, and the protein residues
are in the same order in both files, so :func:`construct_residue_map` transfers the numbering by
position and refuses to guess if the two disagree. Analysing by raw prmtop resid would silently
report the wrong residues.

**Selections.** An Amber prmtop carries no chains and (here) no elements, and Lipid21 has no POPC
residue -- it is a PC headgroup plus PA/OL acyl residues. Everything structural is therefore
selected by mass or by `protein`, never by chain, element or lipid resname. This is the same
lesson, and the same constants, as `02.10.00_slurm_bundle/check_equilibration.py`.

**Interaction geometry.** Contact/H-bond/ionic cutoffs are imported from
:mod:`zh853mor.interactions`, so a contact occupancy measured over MD means the same thing as a
contact in the static cryo-EM fingerprint (Objective 1) and the two can be compared directly.
"""

from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import capped_distance, distance_array

from zh853mor.interactions import CONTACT_CUT, HBOND_CUT
from zh853mor.structure import POCKET_BW

# --- selection constants (masses, not names -- see the module docstring) -----------------------
MASS_P = (30.5, 31.5)     # phosphorus: one per phospholipid, nothing else in this system has it
MASS_O = (15.5, 16.5)     # oxygen, incl. water O (OPC's EPW extra point is massless)
MASS_N = (13.5, 14.5)
MASS_NA = (22.5, 23.5)    # Na+
WATER_RESN = {"WAT", "HOH", "OPC", "TIP3", "SOL"}
STEROL_RESN = {"CHL", "CHL1", "CLR", "CHOL", "CHO"}

CORE_HALF_A = 8.0    # |z - midplane| defining the hydrophobic core
TM_HALF_A = 15.0     # |z - midplane| defining the membrane-embedded CA set
REPORT_PS = 100.0    # DCD cadence written by 03_production.py; fallback if the DCD has no dt

# --- residues followed frame by frame ----------------------------------------------------------
# The ZH853 contact shell (POCKET_BW, from the 4.5 A cryo-EM shell) plus the functional positions
# the pocket shell does not contain. Occupancy is computed for EVERY receptor residue; these are
# the ones whose per-frame distance series is also kept, for equilibration/convergence analysis.
FUNCTIONAL: dict[int, str] = {
    116: "2.50",   # D2.50, the Na+ site -- the ASP/ASH variant under test (SPECIFICATION D-11)
    166: "3.49",   # DRY motif
    167: "3.50",
    168: "3.51",
    281: "6.34",   # TM6 cytoplasmic end; R3.50-T6.34 is the activation ruler
    334: "7.49",   # NPxxY
    335: "7.50",
    338: "7.53",
}
KEY_RESIDUES: dict[int, str] = {**POCKET_BW, **FUNCTIONAL}

# Anchors quoted throughout the project (docs/PLAN.md 1.2): the shared opioid anchors and the two
# candidate ZH853-distinctive contacts. Reported individually, with replicate error bars.
ANCHORS: dict[int, str] = {149: "D3.32", 328: "Y7.43", 299: "H6.52",
                           231: "E-ECL2", 321: "H7.36", 129: "N2.63"}

# Activation-state rulers, as (label, resid_a, resid_b) CA-CA distances. CA rather than sidechain
# tips: at 3.5 A the sidechain rotamers are the least reliable coordinates in the model, and the
# TM6 opening is a backbone motion.
ACTIVATION_DISTANCES: list[tuple[str, int, int]] = [
    ("tm3_tm6_ca", 167, 281),    # R3.50-T6.34: opens ~10 A on activation in class-A GPCRs
    ("npxxy_tm3_ca", 167, 338),  # R3.50-Y7.53: the TM7 half of the activation switch
]

C142_C219 = (142, 219)  # the conserved ECL2-TM3 disulfide, as an integrity check

# tleap renames a residue whose FORM differs from the default without changing WHICH amino acid it
# is: a disulfide cysteine becomes CYX, and a titratable residue carries its protonation state in
# the name. Every build therefore has CYX where the staged receptor.pdb has CYS, at both ends of
# the C142-C219 disulfide. That is not a numbering difference, so the map below folds these to the
# parent residue before comparing sequences -- while a genuinely different amino acid still stops
# the analysis, which is the point of the check.
RESIDUE_FORMS: dict[str, str] = {
    "CYX": "CYS", "CYM": "CYS",
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "TYM": "TYR", "ARN": "ARG",
}


def canonical_residue(name: str) -> str:
    """Residue name with its protonation/disulfide form removed (CYX -> CYS, HID -> HIS)."""
    clean = name.strip().upper()
    return RESIDUE_FORMS.get(clean, clean)


@dataclass(frozen=True)
class BuildDir:
    """One built system: `intermediate/02.10.00_build/<ligand>_<d250>_<timestamp>/`."""

    path: Path
    ligand: str
    d250: str
    stamp: str

    @property
    def name(self) -> str:
        return f"{self.ligand}_{self.d250}"

    @property
    def prmtop(self) -> Path:
        return self.path / "system.prmtop"

    @property
    def receptor_pdb(self) -> Path:
        return self.path / "receptor.pdb"

    def replicas(self) -> list[tuple[str, Path]]:
        """[(replica name, DCD path)] for every production trajectory present, in order."""
        out = []
        for dcd in sorted(self.path.glob("prod_r*.dcd")):
            out.append((dcd.stem, dcd))
        return out

    def sampling(self) -> dict[str, float]:
        """`ZH_*` numbers from this build's `sampling.env`, empty if it has none.

        01_build_system.sh writes the file with the values already substituted (the heredoc is
        unquoted), so a bare `ZH_PROD_NS=500` is the normal case; the `${VAR:-default}` form is
        still accepted in case one was hand-edited back in.
        """
        out: dict[str, float] = {}
        try:
            text = (self.path / "sampling.env").read_text()
        except OSError:
            return out
        for line in text.splitlines():
            m = re.match(r"\s*(ZH_[A-Z_]+)=\$?\{?[A-Z_]*:?-?([0-9]+(?:\.[0-9]+)?)\}?", line)
            if m:
                out[m[1]] = float(m[2])
        return out

    def ligand_resname(self) -> str | None:
        """Resname of the ligand in the built system, or None for an apo build."""
        try:
            meta = json.loads((self.path / "system.json").read_text())
        except (OSError, ValueError):
            return None if self.ligand == "apo" else "LIG"
        name = meta.get("ligand_resname") or ""
        return str(name) or None


_BUILD_RE = re.compile(r"^(?P<ligand>[A-Za-z0-9]+)_(?P<d250>ASP|ASH)_(?P<stamp>\d{8}_\d{6})$")


def discover_builds(root: Path, newest_only: bool = True) -> list[BuildDir]:
    """Build directories under `root`, sorted by (ligand, D2.50 state, time).

    With `newest_only` (the default) a rebuilt system contributes only its most recent directory:
    build directories are timestamped and never reused, so an aborted or superseded build sits
    beside the good one and would otherwise be averaged into it.
    """
    found: list[BuildDir] = []
    for path in sorted(p for p in root.glob("*_*_*") if p.is_dir()):
        m = _BUILD_RE.match(path.name)
        if m and (path / "system.prmtop").exists():
            found.append(BuildDir(path, m["ligand"], m["d250"], m["stamp"]))
    if newest_only:
        latest: dict[str, BuildDir] = {}
        for b in found:
            if b.name not in latest or b.stamp > latest[b.name].stamp:
                latest[b.name] = b
        found = list(latest.values())
    return sorted(found, key=lambda b: (b.ligand != "apo", b.ligand, b.d250, b.stamp))


# --- topology ----------------------------------------------------------------------------------
def protein_residues(u: mda.Universe) -> mda.ResidueGroup:
    """The receptor as an MDAnalysis ResidueGroup, in topology order."""
    prot = u.select_atoms("protein")
    if not len(prot):
        raise ValueError("no atoms matched `protein`; is this a receptor system?")
    return prot.residues


def construct_residue_map(u: mda.Universe, receptor_pdb: Path,
                          notes: list[str] | None = None) -> np.ndarray:
    """Human-OPRM1 residue id for each `protein` residue of `u`, in topology order.

    tleap renumbers from 1; `receptor.pdb` (the staged, OPM-oriented input the system was built
    from) keeps the construct numbering. The two contain the same protein residues in the same
    order, so the map is positional. Verified by residue name, not assumed: a mismatch means the
    build directory's receptor.pdb is not the file the prmtop was built from, and every residue
    number downstream would be wrong in a way that still produces plausible numbers.

    Names are compared after folding away protonation and disulfide FORM (see `RESIDUE_FORMS`),
    which tleap assigns and the staged PDB does not carry: those differ in every build and say
    nothing about the numbering. Any that are found are appended to `notes` -- expected for the
    C142-C219 cysteines, worth reading for a histidine, since a tautomer that disagrees with the
    prepared receptor is the D-15 failure returning.
    """
    sys_res = protein_residues(u)
    ref_res = protein_residues(mda.Universe(str(receptor_pdb)))
    sys_names = [r.resname.strip().upper() for r in sys_res]
    ref_names = [r.resname.strip().upper() for r in ref_res]
    if len(sys_names) != len(ref_names):
        raise ValueError(
            f"{receptor_pdb.name} has {len(ref_names)} protein residues but the topology has "
            f"{len(sys_names)}; they are not the same receptor."
        )
    resids = np.array([int(r.resid) for r in ref_res], dtype=int)
    bad, forms = [], []
    for i, (a, b) in enumerate(zip(sys_names, ref_names, strict=True)):
        if a == b or canonical_residue(a) == canonical_residue(b) == "CYS":
            # A disulfide cysteine is CYX in every prmtop this pipeline builds and CYS in every
            # staged receptor. Reporting it each run would be noise around the cases that matter.
            continue
        if canonical_residue(a) == canonical_residue(b):
            forms.append(f"{b}{resids[i]}->{a}")
        else:
            bad.append((int(resids[i]), a, b))
    if bad:
        resid, a, b = bad[0]
        raise ValueError(
            f"residue {resid} is {a} in the topology but {b} in {receptor_pdb.name} "
            f"({len(bad)} mismatches); the numbering transfer would be wrong."
        )
    if forms and notes is not None:
        # Not the disulfides: a PROTONATION form that disagrees with the prepared receptor is the
        # D-15 failure (tleap overriding an upstream tautomer) coming back, and is worth reading.
        notes.append(f"{len(forms)} residue(s) carry a different protonation form in the topology "
                     f"than in {receptor_pdb.name}: {', '.join(forms)} -- check that tleap did not "
                     "override the prepared assignment (SPECIFICATION D-15)")
    return resids


def construct_ids_for(atoms: mda.AtomGroup, prot: mda.AtomGroup,
                      resids: np.ndarray) -> np.ndarray:
    """Construct residue id for each residue represented in `atoms`, in topology order.

    `resids` is the `construct_residue_map` array, one entry per residue of `prot`. A selection
    is rarely one-atom-per-residue -- the ACE/NME caps carry no CA, and a glycine no sidechain
    polar atom -- so any array derived from a selection needs its OWN residue axis. Indexing such
    an array with the full `resids` silently mis-labels every entry after the first gap (or, if
    you are lucky, raises a shape error in a plot).
    """
    position = {int(r.resid): i for i, r in enumerate(prot.residues)}
    return np.array([resids[position[int(r.resid)]] for r in atoms.residues], dtype=int)


def select_by_mass(u: mda.Universe, lo: float, hi: float, extra: str = "") -> mda.AtomGroup:
    """Atoms whose mass lies in (lo, hi) -- the portable way to select an element here."""
    sel = f"prop mass > {lo} and prop mass < {hi}"
    return u.select_atoms(f"({sel}) and ({extra})" if extra else sel)


def count_lipids(u: mda.Universe) -> tuple[int, int]:
    """(phospholipids, sterols) as MOLECULES, not residues.

    Lipid21 splits one POPC into three residues (PC + PA + OL), so counting residues overstates
    the lipid count threefold and the area per lipid comes out three times too small.
    """
    n_phos = len(select_by_mass(u, *MASS_P))
    n_sterol = sum(1 for r in u.residues if r.resname.strip().upper() in STEROL_RESN)
    return n_phos, n_sterol


# --- per-frame geometry --------------------------------------------------------------------
class ResidueMinDistance:
    """Per-residue minimum distance from a fixed atom set to a moving partner group.

    Built once and called per frame. The whole distance matrix is computed in one C call and
    reduced per residue with `reduceat`, which is what makes a 5000-frame trajectory tractable;
    looping over residues in Python is ~300x slower for the same numbers.
    """

    def __init__(self, atoms: mda.AtomGroup) -> None:
        self.atoms = atoms
        resindices = np.asarray(atoms.resindices)
        if resindices.size and np.any(np.diff(resindices) < 0):
            raise ValueError("atoms must be in topology order for per-residue reduction")
        self.offsets = np.flatnonzero(np.r_[True, np.diff(resindices) != 0])
        self.resindices = resindices[self.offsets]

    def __call__(self, partner_positions: np.ndarray, box: np.ndarray | None) -> np.ndarray:
        if not len(self.atoms) or not len(partner_positions):
            return np.full(self.offsets.size, np.inf)
        d = distance_array(self.atoms.positions, np.asarray(partner_positions, dtype=np.float32),
                           box=box)
        return np.minimum.reduceat(d.min(axis=1), self.offsets)


def occupancy(distances: np.ndarray, cutoff: float = CONTACT_CUT) -> np.ndarray:
    """Fraction of frames within `cutoff`, for a (n_frames, n_residues) distance matrix."""
    return (np.asarray(distances) <= cutoff).mean(axis=0)


def water_bridges(water_o_positions: np.ndarray, residue_polar: list[np.ndarray],
                  ligand_polar: np.ndarray, box: np.ndarray | None,
                  cutoff: float = HBOND_CUT) -> np.ndarray:
    """Waters simultaneously within `cutoff` of the ligand and of each listed residue.

    One count per entry of `residue_polar`. The canonical MOR H6.52 phenol contact is
    water-mediated (docs/PLAN.md 1.2 puts H299 at 4.89 A, outside the direct shell), so a
    direct-contact occupancy alone would score it as absent -- this is the measurement that
    distinguishes "not interacting" from "interacting through a bridging water".
    """
    n = len(residue_polar)
    if not len(water_o_positions) or not len(ligand_polar):
        return np.zeros(n, dtype=int)
    pairs = capped_distance(np.asarray(ligand_polar, dtype=np.float32),
                            np.asarray(water_o_positions, dtype=np.float32),
                            max_cutoff=cutoff, box=box, return_distances=False)
    near_lig = np.unique(np.asarray(pairs)[:, 1]) if len(pairs) else np.array([], dtype=int)
    if not near_lig.size:
        return np.zeros(n, dtype=int)
    shortlist = np.asarray(water_o_positions, dtype=np.float32)[near_lig]
    out = np.zeros(n, dtype=int)
    for i, polar in enumerate(residue_polar):
        if polar is None or not len(polar):
            continue
        d = distance_array(np.asarray(polar, dtype=np.float32), shortlist, box=box)
        out[i] = int((d.min(axis=0) <= cutoff).sum())
    return out


def membrane_midplane(phosphorus_positions: np.ndarray) -> float:
    """Bilayer midplane z as the mean phosphate height."""
    return float(np.mean(np.asarray(phosphorus_positions)[:, 2]))


def bilayer_thickness(phosphorus_positions: np.ndarray) -> float:
    """Phosphate-to-phosphate thickness (A): upper-leaflet mean minus lower-leaflet mean."""
    z = np.asarray(phosphorus_positions)[:, 2]
    if z.size < 2:
        return float("nan")
    mid = float(z.mean())
    up, lo = z[z >= mid], z[z < mid]
    if not up.size or not lo.size:
        return float("nan")
    return float(up.mean() - lo.mean())


def protein_cross_section(positions: np.ndarray, midplane_z: float,
                          half_width: float = TM_HALF_A) -> float | None:
    """xy convex-hull area (A^2) of the membrane-embedded protein, or None without SciPy.

    A hull overestimates a non-convex cross-section, so an area per lipid derived from it is a
    LOWER bound -- which is why both the gross and the corrected values are reported rather than
    one replacing the other. (Same definition as `check_equilibration.py`, so the production
    numbers can be read against the equilibration report.)
    """
    try:
        from scipy.spatial import ConvexHull
    except ImportError:
        return None
    sl = np.asarray(positions)[np.abs(np.asarray(positions)[:, 2] - midplane_z) < half_width][:, :2]
    if len(sl) < 3:
        return None
    return float(ConvexHull(sl).volume)  # for 2-D input, `volume` is the area


def area_per_lipid(box_xy: tuple[float, float], n_lipid_molecules: int,
                   protein_area: float = 0.0) -> float:
    """Area per lipid (A^2): (box cross-section - protein cross-section) / lipids per leaflet.

    With `protein_area` left at 0 this is the GROSS value, systematically high for a
    protein-containing bilayer: it is then a drift observable rather than a number to compare
    against a pure-POPC literature value.
    """
    per_leaflet = max(n_lipid_molecules / 2.0, 1.0)
    return float((box_xy[0] * box_xy[1] - protein_area) / per_leaflet)


def principal_components(coords: np.ndarray, n_components: int = 3) -> tuple[np.ndarray, np.ndarray]:
    """(projections (n_frames, k), explained-variance ratios) of aligned Cartesian coordinates.

    `coords` is (n_frames, n_atoms, 3) and must already be superposed -- otherwise the leading
    components are rigid-body motion. Feeds `convergence.cosine_content`.
    """
    x = np.asarray(coords, dtype=float)
    flat = x.reshape(x.shape[0], -1)
    flat = flat - flat.mean(axis=0)
    # SVD of the (frames x 3N) matrix rather than an explicit 3N x 3N covariance: for 281 CA
    # (843 columns) and 5000 frames this is both faster and better conditioned.
    u_, s, _ = np.linalg.svd(flat, full_matrices=False)
    k = min(n_components, s.size)
    var = s**2
    return u_[:, :k] * s[:k], (var[:k] / var.sum()) if var.sum() > 0 else np.zeros(k)


@dataclass
class ReplicaSpan:
    """How much trajectory a replica actually contains, and whether that is all of it.

    Read from the DCD HEADER and the state log, never by loading frames, so listing a whole panel
    costs milliseconds: an unfinished or dead replica should be visible before a job array is
    submitted, not after it has been reduced.
    """

    n_frames: int
    dt_ps: float
    ns: float
    target_ns: float | None
    log_ns: float | None
    age_s: float          # seconds since the DCD was last written

    @property
    def fraction(self) -> float | None:
        if not self.target_ns:
            return None
        return self.ns / self.target_ns

    @property
    def status(self) -> str:
        """One word for the listing: what state this replica is in."""
        if self.age_s < WRITING_WINDOW_S:
            return "writing"
        if self.target_ns is None:
            return "no-target"
        frac = self.fraction or 0.0
        if frac >= COMPLETE_FRACTION:
            return "complete"
        return f"partial-{frac * 100:.0f}%"


WRITING_WINDOW_S = 900.0    # a DCD touched within 15 min is probably still being written
COMPLETE_FRACTION = 0.99    # the last frame lands one report interval short of the nominal length


def dcd_span(dcd: Path) -> tuple[int, float]:
    """(frames, ps per frame) from the DCD header alone -- no topology, no frames read.

    The DCD header carries the frame count and the inter-frame interval (AKMA), so this is O(1)
    on a 5 GB trajectory.
    """
    from MDAnalysis import units
    from MDAnalysis.lib.formats.libdcd import DCDFile

    with DCDFile(str(dcd)) as fh:
        header = dict(fh.header)
        n_frames = int(fh.n_frames)
    dt_ps = float(units.convert(float(header["delta"]), "AKMA", "ps"))
    dt_ps *= max(int(header.get("nsavc", 1)), 1)
    return n_frames, dt_ps


def replica_span(dcd: Path, target_ns: float | None = None) -> ReplicaSpan:
    """How long the replica at `dcd` ran, against the length its build was configured for."""
    n_frames, dt_ps = dcd_span(dcd)
    log_ns = None
    log = dcd.with_suffix(".log")
    if log.exists():
        try:
            series = read_state_log(log)
            if "time_ps" in series and series["time_ps"].size:
                log_ns = float(series["time_ps"][-1]) / 1000.0
        except (OSError, ValueError):
            log_ns = None
    return ReplicaSpan(n_frames=n_frames, dt_ps=dt_ps, ns=n_frames * dt_ps / 1000.0,
                       target_ns=target_ns, log_ns=log_ns,
                       age_s=max(time.time() - dcd.stat().st_mtime, 0.0))


def read_state_log(path: Path) -> dict[str, np.ndarray]:
    """Parse an OpenMM StateDataReporter CSV BY COLUMN NAME.

    The reporter emits its columns in a fixed internal order, not the order the kwargs were
    passed, so positional parsing silently mislabels energy as temperature the moment a field is
    added. (Same parser as `check_equilibration.py`; the thermodynamic series are the half of
    production QC that needs no trajectory.)
    """
    import csv

    with open(path) as fh:
        rows = list(csv.reader(fh))
    if len(rows) < 2:
        raise ValueError(f"{path} has no data rows; did the run die early?")
    header = [h.lstrip("#").strip().strip('"') for h in rows[0]]
    data = np.array([[float(x) for x in r] for r in rows[1:] if len(r) == len(header)])
    if data.size == 0:
        raise ValueError(f"{path} has a header but no parsable rows.")

    def col(*needles: str) -> np.ndarray | None:
        for i, h in enumerate(header):
            if all(n.lower() in h.lower() for n in needles):
                return data[:, i]
        return None

    wanted = {"step": ("Step",), "time_ps": ("Time",), "pe": ("Potential", "Energy"),
              "temperature": ("Temperature",), "density": ("Density",),
              "volume": ("Box", "Volume"), "speed": ("Speed",)}
    out = {}
    for key, needles in wanted.items():
        found = col(*needles)
        if found is not None:
            out[key] = found
    return out


@dataclass
class ReplicaResult:
    """One reduced replica: the JSON summary plus lazy access to its timeseries arrays."""

    summary: dict
    npz_path: Path

    @property
    def system(self) -> str:
        return str(self.summary.get("build", self.npz_path.parent.name))

    @property
    def replica(self) -> str:
        return str(self.summary.get("replica", self.npz_path.stem))

    @property
    def ligand(self) -> str:
        return str(self.summary.get("ligand", self.system.split("_")[0]))

    @property
    def d250(self) -> str:
        return str(self.summary.get("d250", self.system.rsplit("_", 1)[-1]))

    @property
    def t0(self) -> int:
        eq = self.summary.get("equilibration", {})
        return int(eq.get("t0_frames", 0))

    def arrays(self) -> dict[str, np.ndarray]:
        """Load the .npz. Kept explicit rather than cached: the caller decides what to hold."""
        with np.load(self.npz_path, allow_pickle=False) as data:
            return {k: data[k] for k in data.files}

    def series(self, name: str, equilibrated: bool = True) -> np.ndarray:
        """One timeseries, by default with the pre-equilibration frames dropped."""
        with np.load(self.npz_path, allow_pickle=False) as data:
            if name not in data.files:
                raise KeyError(f"{self.npz_path.name} has no series '{name}'")
            arr = np.asarray(data[name], dtype=float)
        return arr[self.t0:] if equilibrated else arr

    def mean_of(self, name: str) -> float:
        """Post-equilibration mean of an observable, as recorded by the reduction step."""
        entry = self.summary.get("means", {}).get(name)
        return float(entry["mean"]) if entry else float("nan")


def load_replicas(root: Path) -> list[ReplicaResult]:
    """Every reduced replica under `root` (intermediate/02.11.00_analyze_simulations/<system>/<rep>.json)."""
    out = []
    for js in sorted(root.glob("*/*.json")):
        npz = js.with_suffix(".npz")
        if not npz.exists():
            continue
        try:
            summary = json.loads(js.read_text())
        except ValueError:
            continue
        out.append(ReplicaResult(summary, npz))
    return out


def group_by_system(results: list[ReplicaResult]) -> dict[str, list[ReplicaResult]]:
    """Replicas keyed by system name, apo first and then ligands alphabetically.

    Ordering is fixed here so every table and figure in the stage lists the systems the same way.
    """
    systems: dict[str, list[ReplicaResult]] = {}
    for r in results:
        systems.setdefault(r.system, []).append(r)
    order = sorted(systems, key=lambda s: (not s.startswith("apo"), s))
    return {s: sorted(systems[s], key=lambda r: r.replica) for s in order}

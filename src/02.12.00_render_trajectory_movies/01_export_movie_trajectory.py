#!/usr/bin/env python3
"""Export one production replica as a small, viewable movie trajectory.

**Runs on the cluster** (`zh853mor-prep`, CPU), beside `02.11.00`'s reduction and for the same
reason: a 500 ns replica is ~5 GB and the panel is ~100 GB, so the trajectory never leaves the
machine it was written on. This step writes the few MB per replica that a movie actually needs --
a topology PDB, an XTC, and a JSON of provenance -- which is what gets copied back.

What "movie-ready" means here, and why each part is necessary:

* **Whole molecules.** OpenMM writes coordinates wrapped into the periodic box molecule by
  molecule. Rendered as-is, lipids and the receptor are torn in half at the box faces and the
  ligand can appear on the far side of the box while still sitting in the pocket. Every frame is
  therefore unwrapped by fragment, re-centred on the protein, and re-wrapped around it.
* **A fixed atom selection.** A trajectory file has one atom count for every frame, so a
  "lipids within 8 A" selection cannot be re-evaluated per frame. The annulus is chosen ONCE,
  as the union over frames sampled across the whole run, and then held fixed.
* **Superposition.** Rotational and translational diffusion of the whole system is real but
  uninformative, and it dominates a movie. Frames are superposed on the membrane-embedded CA,
  the same set `02.11.00` measures its TM RMSD on. `--no-align` keeps the raw box view, which is
  what you want when the question is tilt, drift or a PBC artefact rather than internal motion.

The output pair opens directly in VMD, PyMOL, ChimeraX and MolStar, and is what
`02_render_dashboard.py` (matplotlib diagnostics) and `03_render_molstar.js` (cartoon render)
both read.

    python 01_export_movie_trajectory.py --list
    python 01_export_movie_trajectory.py --build ../../intermediate/02.10.00_build/apo_ASH_20260901_133002
    ./submit_export.sh                      # one array task per replica
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

import MDAnalysis as mda  # noqa: E402
import numpy as np  # noqa: E402
from MDAnalysis import transformations as trans  # noqa: E402
from MDAnalysis.analysis.align import rotation_matrix  # noqa: E402
from MDAnalysis.lib.distances import capped_distance  # noqa: E402

from zh853mor import md, paths  # noqa: E402

BUILD_ROOT = paths.INTERMEDIATE / "02.10.00_build"
OUT_ROOT = paths.INTERMEDIATE / "02.12.00_render_trajectory_movies"

# --- what goes into the movie ------------------------------------------------------------------
# Defaults chosen so a replica lands in the low tens of MB: small enough to scp back and to load
# in a browser, large enough to show the things a time series cannot.
MOVIE_FRAMES = 300       # ~12 s at 25 fps; the eye cannot follow more in one pass anyway
SCAN_FRAMES = 40         # frames sampled to CHOOSE the fixed selection (see module docstring)
LIPID_CUT = 8.0          # A; annular shell -- the lipids that actually touch the receptor
WATER_CUT = 5.0          # A; a water this close to the ligand can bridge a polar contact
N_WATERS = 40            # the most persistent pocket waters, by occupancy over the scan frames
MAX_ATOMS = 40_000       # guard: past this an XTC frame stops being something a browser enjoys

# Group order in the written files. Each group is a CONTIGUOUS slice of the movie atom index, so
# the JSON can describe it as [start, stop) and the renderers need no selection language.
GROUP_ORDER = ("protein", "ligand", "lipid", "phosphate", "water", "ion")

# Chain IDs written into the movie PDB. The non-protein ones are fixed so a renderer can select
# "the ligand" as `chain L` without parsing residue names, and the protein alphabet deliberately
# EXCLUDES them -- a fifth protein chain landing on L would otherwise be drawn as the ligand.
CHAIN_LIGAND, CHAIN_LIPID, CHAIN_WATER, CHAIN_ION = "L", "M", "W", "I"
PROTEIN_CHAINS = "".join(c for c in "ABCDEFGHJKNOPQRSTUVXYZ"
                         if c not in (CHAIN_LIGAND, CHAIN_LIPID, CHAIN_WATER, CHAIN_ION))


@dataclass(frozen=True)
class Job:
    build: md.BuildDir
    replica: str
    dcd: Path


def jobs(builds: list[md.BuildDir]) -> list[Job]:
    return [Job(b, name, dcd) for b in builds for name, dcd in b.replicas()]


def listing(todo: list[Job], frames: int) -> list[str]:
    """One tab-separated row per job: what it is, how long it is, and what it will produce.

    Same shape as `02.11.00`'s listing -- and for the same reason: a partial run, a run still
    being written and a finished one are indistinguishable by file name, and exporting the first
    two silently makes a movie of a run that is not what it claims to be. The extra columns are
    the ones specific to this stage: how many frames survive the stride, and at what spacing in
    simulated time, because a 500 ns replica squeezed into 300 frames is sampled every 1.7 ns and
    anything faster than that is invisible in the result.
    """
    rows = []
    for i, job in enumerate(todo, start=1):
        target = job.build.sampling().get("ZH_PROD_NS")
        span = md.replica_span(job.dcd, target_ns=target)
        stride = max(1, span.n_frames // frames) if span.n_frames else 1
        n_out = len(range(0, span.n_frames, stride)[:frames]) if span.n_frames else 0
        per_frame = stride * span.dt_ps / 1000.0
        rows.append("\t".join([
            str(i), job.build.name, job.replica,
            str(span.n_frames) if span.n_frames else "?",
            f"{span.ns:.1f}" if span.n_frames else "?",
            str(n_out), f"{per_frame:.2f}" if n_out else "?", f"1/{stride}",
            span.status,
        ]))
    return rows


def inventory_notes(builds: list[md.BuildDir]) -> list[str]:
    """Remarks about whole builds rather than single replicas (chiefly: a missing replica)."""
    notes = []
    for b in builds:
        present = {name for name, _ in b.replicas()}
        want = int(b.sampling().get("ZH_REPLICAS", 0))
        if want and len(present) < want:
            notes.append(f"note: {b.name} has {len(present)} of {want} configured replicas")
    return notes


# --- topology grouping --------------------------------------------------------------------------
def classify(u: mda.Universe, ligand_resname: str | None,
             warnings: list[str]) -> dict[str, mda.AtomGroup]:
    """Split the system into protein / ligand / lipid / water / ion, by composition not by name.

    The only name this depends on is the ligand's, which `system.json` records. Everything else is
    identified structurally, the way the rest of the package does it (`zh853mor.md`): water by
    residue name from a known set, ions as the MONATOMIC fragments left over, lipids as whatever
    else is not protein. Lipid21 splits a POPC across three residues (head + two tails) bonded
    into one molecule, so anything lipid-related must work on fragments, never on residues.
    """
    protein = u.select_atoms("protein")
    if not len(protein):
        raise ValueError("no atoms matched `protein`; is this a receptor system?")

    water = u.select_atoms("resname " + " ".join(sorted(md.WATER_RESN)))
    ligand = u.select_atoms(f"resname {ligand_resname}") if ligand_resname else u.atoms[[]]
    if ligand_resname and not len(ligand):
        warnings.append(f"system.json declares ligand resname {ligand_resname} but no atoms "
                        "match it; the movie is of an apo-looking system")

    rest = u.atoms - protein - water - ligand
    # An ion is a molecule of one atom. Defining it that way rather than by mass or name catches
    # the Cl- counter-ions as well as the Na+ the D2.50 analysis is about, and cannot mistake a
    # lipid phosphorus for one.
    ion_ix, lipid_ix = [], []
    for frag in rest.fragments:
        (ion_ix if len(frag) == 1 else lipid_ix).extend(frag.ix)
    return {
        "protein": protein,
        "ligand": ligand,
        "lipid": u.atoms[np.array(sorted(lipid_ix), dtype=int)] if lipid_ix else u.atoms[[]],
        "water": water,
        "ion": u.atoms[np.array(sorted(ion_ix), dtype=int)] if ion_ix else u.atoms[[]],
    }


def centre_laterally(protein: mda.AtomGroup, centre_z: bool):
    """Put the protein at the middle of the box IN X AND Y, and in z only if asked.

    x and y are centred unconditionally: the box origin is arbitrary in the membrane plane, and
    without this the receptor wanders out of frame and its annular lipids are imaged to the far
    face. z is left alone unless the frames are going to be superposed anyway, because along the
    membrane normal the position IS meaningful -- a receptor sliding out of the bilayer, or a
    bilayer drifting through the box, is one of the things an unaligned movie is watched for.
    """
    def transform(ts):
        cog = protein.positions.mean(axis=0)
        shift = np.zeros(3, dtype=np.float32)
        shift[:2] = np.asarray(ts.dimensions[:2], dtype=float) / 2.0 - cog[:2]
        if centre_z:
            shift[2] = float(ts.dimensions[2]) / 2.0 - cog[2]
        ts.positions += shift
        return ts
    return transform


def imaged(u: mda.Universe, protein: mda.AtomGroup, centre_z: bool,
           warnings: list[str]) -> None:
    """Install the unwrap/centre/rewrap workflow so EVERY frame read is already whole.

    Attached to the reader rather than applied by hand so that the frames scanned to choose the
    selection and the frames written to the XTC are imaged identically -- if they were not, a
    lipid picked as "annular" in the scan could be written a box-length away.
    """
    if not hasattr(u, "bonds") or len(u.bonds) == 0:
        # A prmtop always has bonds; a PDB passed via --topology for a smoke test may not, and
        # without them `unwrap` cannot tell which atoms belong to the same molecule.
        warnings.append("topology carries no bonds; guessing them by distance (slow, and only "
                        "sound for a small test system)")
        u.atoms.guess_bonds()
    u.trajectory.add_transformations(
        trans.unwrap(u.atoms),
        centre_laterally(protein, centre_z),
        trans.wrap(u.atoms, compound="fragments"),
    )


def tm_reference(u: mda.Universe, protein: mda.AtomGroup,
                 warnings: list[str]) -> tuple[mda.AtomGroup, str]:
    """The CA used for superposition: membrane-embedded only, as in `02.11.00`.

    Superposing on every CA lets a floppy terminus or, in system A, the whole Gi heterotrimer
    pull the receptor around in the frame. The TM bundle is the part that is supposed to hold
    still, so it is the part the camera should be bolted to.
    """
    ca = protein.select_atoms("name CA")
    everything = f"all {len(ca)} CA"
    phos = md.select_by_mass(u, *md.MASS_P)
    if not len(phos):
        warnings.append("no phosphorus found: superposing on every CA rather than the TM subset")
        return ca, everything
    mid = md.membrane_midplane(phos.positions)
    tm = ca[np.abs(ca.positions[:, 2] - mid) < md.TM_HALF_A]
    if len(tm) < 100:
        warnings.append(f"only {len(tm)} CA inside the membrane slab; superposing on every CA")
        return ca, everything
    return tm, f"{len(tm)} membrane-embedded CA (|z - midplane| < {md.TM_HALF_A} A)"


# --- choosing the fixed selection ---------------------------------------------------------------
def scan_selection(u: mda.Universe, groups: dict[str, mda.AtomGroup], frame_ix: np.ndarray,
                   lipid_cut: float, water_cut: float, n_waters: int, pocket: mda.AtomGroup,
                   ) -> tuple[np.ndarray, np.ndarray, list[tuple[int, float]]]:
    """Pick the annular lipids and the persistent pocket waters, over `frame_ix`.

    Returns the lipid ATOM indices (whole fragments, union over the scanned frames) and the
    waters as (resid, occupancy) sorted by occupancy. Lipids are taken as a union because a lipid
    that visits the annulus is worth seeing arrive; waters are ranked and truncated instead,
    because over 500 ns thousands of waters pass within 5 A of the pocket and drawing them all
    would hide the handful that stay long enough to mediate a contact.
    """
    protein_heavy = groups["protein"].select_atoms("prop mass > 2.0")
    lipids, water = groups["lipid"], groups["water"]
    # fragindex per atom, so a hit on any atom pulls in the whole lipid molecule -- necessary
    # because Lipid21 spreads one POPC over three residues.
    lipid_frag = lipids.fragindices if len(lipids) else np.empty(0, int)
    water_res = water.resindices if len(water) else np.empty(0, int)

    hit_frags: set[int] = set()
    water_hits: dict[int, int] = {}
    for fi in frame_ix:
        ts = u.trajectory[int(fi)]
        box = ts.dimensions
        if len(lipids) and lipid_cut > 0:
            pairs = capped_distance(protein_heavy.positions, lipids.positions,
                                    max_cutoff=lipid_cut, box=box, return_distances=False)
            hit_frags.update(np.unique(lipid_frag[pairs[:, 1]]).tolist())
        if len(water) and len(pocket) and n_waters > 0:
            pairs = capped_distance(pocket.positions, water.positions,
                                    max_cutoff=water_cut, box=box, return_distances=False)
            for ri in np.unique(water_res[pairs[:, 1]]):
                water_hits[int(ri)] = water_hits.get(int(ri), 0) + 1

    lipid_ix = (lipids[np.isin(lipid_frag, list(hit_frags))].ix
                if hit_frags else np.empty(0, dtype=int))
    ranked = sorted(water_hits.items(), key=lambda kv: -kv[1])[:n_waters]
    waters = [(int(u.residues[ri].resid), n / len(frame_ix)) for ri, n in ranked]
    water_ix = (np.concatenate([u.residues[ri].atoms.ix for ri, _ in ranked])
                if ranked else np.empty(0, dtype=int))
    return lipid_ix, water_ix, waters


def assemble(u: mda.Universe, groups: dict[str, mda.AtomGroup], lipid_ix: np.ndarray,
             water_ix: np.ndarray, hydrogens: bool, membrane: bool,
             warnings: list[str]) -> tuple[np.ndarray, dict[str, list[int]]]:
    """Concatenate the kept atoms in `GROUP_ORDER`, disjointly, and report each group's slice."""
    def heavy(ix: np.ndarray) -> np.ndarray:
        if hydrogens or not len(ix):
            return ix
        return ix[u.atoms[ix].masses > 2.0]

    # Phosphorus WITHIN the lipids, not wherever a mass of ~31 turns up: with a prmtop the two
    # are the same set, but a topology whose masses were guessed (a PDB read back for a smoke
    # test) puts OPC's EPW extra point at 30.97 and would dot every water as a lipid headgroup.
    phosphate = (np.intersect1d(md.select_by_mass(u, *md.MASS_P).ix, groups["lipid"].ix)
                 if membrane else np.empty(0, dtype=int))
    proposed: dict[str, np.ndarray] = {
        "protein": heavy(groups["protein"].ix),
        "ligand": heavy(groups["ligand"].ix),
        "lipid": heavy(lipid_ix),
        # The rest of the bilayer as phosphorus only: two dotted leaflets are enough to read the
        # membrane plane, undulation and the receptor's tilt against it, for ~300 atoms.
        "phosphate": phosphate,
        # Water hydrogens are kept whatever --hydrogens says: a lone oxygen sphere cannot show
        # which way a bridging water is pointing, which is the only reason these waters are here.
        # OPC's massless EPW extra point is dropped -- it is a charge site, not an atom.
        "water": water_ix[u.atoms[water_ix].masses > 0.0] if len(water_ix) else water_ix,
        "ion": groups["ion"].ix,
    }
    # The groups must be DISJOINT: an atom written twice breaks the [start, stop) slices the
    # renderers colour by, and MDAnalysis would happily write it twice into the XTC. Earlier
    # groups win, in GROUP_ORDER -- the annular lipid keeps its phosphorus, the dotted plane
    # does not repeat it.
    parts: dict[str, np.ndarray] = {}
    taken = np.empty(0, dtype=int)
    for g in GROUP_ORDER:
        parts[g] = np.setdiff1d(proposed[g], taken) if len(taken) else np.asarray(proposed[g])
        taken = np.union1d(taken, parts[g])
    index = np.concatenate([parts[g] for g in GROUP_ORDER])
    if len(index) != len(np.unique(index)):
        raise AssertionError("movie selection contains duplicated atoms")

    slices, at = {}, 0
    for g in GROUP_ORDER:
        slices[g] = [at, at + len(parts[g])]
        at += len(parts[g])
    if len(index) > MAX_ATOMS:
        warnings.append(f"{len(index)} atoms is above the {MAX_ATOMS} guard; the XTC and the "
                        "browser render will both be slow. Lower --lipid-cut or --waters.")
    return index, slices


# --- the export ---------------------------------------------------------------------------------
def export_replica(job: Job, out_dir: Path, *, frames: int = MOVIE_FRAMES,
                   scan_frames: int = SCAN_FRAMES, lipid_cut: float = LIPID_CUT,
                   water_cut: float = WATER_CUT, n_waters: int = N_WATERS,
                   hydrogens: bool = False, membrane: bool = True, align: bool = True,
                   topology: Path | None = None) -> dict[str, object]:
    """Write `<replica>_movie.{pdb,xtc,json}` for one replica and return the JSON summary."""
    started = time.time()
    build = job.build
    warnings: list[str] = []
    top = topology or build.prmtop
    u = mda.Universe(str(top), str(job.dcd))

    resname = build.ligand_resname()
    groups = classify(u, resname, warnings)
    imaged(u, groups["protein"], centre_z=align, warnings=warnings)

    n_src = len(u.trajectory)
    stride = max(1, n_src // frames)
    keep = np.arange(0, n_src, stride)[:frames]
    dt_ps = float(u.trajectory.dt) or md.REPORT_PS
    if not np.isfinite(dt_ps) or dt_ps <= 0:
        dt_ps = md.REPORT_PS
        warnings.append(f"trajectory reports no timestep; assuming {md.REPORT_PS} ps per frame")
    time_ns = keep * dt_ps / 1000.0

    # The pocket, for the water scan: the ligand if there is one, otherwise the receptor residues
    # that would line it. An apo run's pocket waters are exactly what the holo comparison needs.
    resids, numbering = construct_numbering(u, groups["protein"], build, warnings)
    if len(groups["ligand"]):
        pocket = groups["ligand"]
    else:
        wanted = [i for i, r in enumerate(resids) if int(r) in md.POCKET_BW]
        pocket = (groups["protein"].residues[wanted].atoms.select_atoms("prop mass > 2.0")
                  if wanted else u.atoms[[]])
        if not len(pocket):
            warnings.append("no ligand and no pocket residues identified; no waters are kept")

    # Sampled ACROSS the run, not from its start: a lipid that only reaches the annulus at 400 ns
    # is exactly the event a movie is for, and a scan of the first frames would never see it.
    scan_ix = keep[np.unique(np.linspace(0, len(keep) - 1, min(scan_frames, len(keep))).astype(int))]
    lipid_ix, water_ix, waters = scan_selection(
        u, groups, scan_ix, lipid_cut, water_cut, n_waters, pocket)
    index, slices = assemble(u, groups, lipid_ix, water_ix, hydrogens, membrane, warnings)
    subset = u.atoms[index]

    # Back to the first movie frame before fixing the TM subset and the superposition reference:
    # the scan left the reader on its last sampled frame, and the membrane slab that defines the
    # TM set must be read at the frame the rest of the movie is fitted to.
    u.trajectory[int(keep[0])]
    tm, align_label = tm_reference(u, groups["protein"], warnings)
    ref = tm.positions.copy()
    ref_com = ref.mean(axis=0)
    ref -= ref_com

    out_dir = paths.ensure_dir(out_dir)
    stem = f"{job.replica}_movie"
    xtc_path, pdb_path = out_dir / f"{stem}.xtc", out_dir / f"{stem}.pdb"
    boxes, fit_rmsd, chains = [], [], {}
    with mda.Writer(str(xtc_path), n_atoms=len(subset)) as W:
        for n, fi in enumerate(keep):
            ts = u.trajectory[int(fi)]
            if align:
                mob = tm.positions
                com = mob.mean(axis=0)
                R, rms = rotation_matrix(mob - com, ref)
                # Rotate about the fitted centre and put it back where the reference sat, so the
                # receptor stays put in the frame instead of sliding with its own centroid.
                subset.positions = (subset.positions - com) @ R.T + ref_com
                fit_rmsd.append(float(rms))
            boxes.append([float(x) for x in ts.dimensions[:3]])
            if n == 0:
                chains = write_topology(pdb_path, subset, resids, groups, slices, numbering)
            W.write(subset)

    summary: dict[str, object] = {
        "build": build.name, "replica": job.replica,
        "ligand": build.ligand, "d250": build.d250,
        "source": {"dcd": str(job.dcd), "topology": str(top), "n_frames": int(n_src),
                   "dt_ps": dt_ps, "length_ns": float(n_src * dt_ps / 1000.0)},
        "movie": {"n_frames": int(len(keep)), "stride": int(stride),
                  "frame_spacing_ns": float(stride * dt_ps / 1000.0),
                  "time_ns": [round(float(t), 4) for t in time_ns],
                  "aligned": bool(align),
                  "align_selection": align_label if align else "none",
                  "fit_rmsd_a": [round(r, 3) for r in fit_rmsd],
                  "box_xyz": boxes},
        "atoms": {"total": int(len(subset)), "groups": slices, "chains": chains,
                  "hydrogens": bool(hydrogens), "numbering": numbering},
        "selection": {"lipid_cut_a": lipid_cut, "water_cut_a": water_cut,
                      "n_waters_kept": len(waters), "scan_frames": int(len(scan_ix)),
                      "membrane": bool(membrane), "ligand_resname": resname},
        "waters": [{"resid": r, "occupancy": round(o, 3)} for r, o in waters],
        "files": {"topology": pdb_path.name, "trajectory": xtc_path.name},
        "warnings": warnings,
        "wall_seconds": round(time.time() - started, 1),
    }
    (out_dir / f"{stem}.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def construct_numbering(u: mda.Universe, protein: mda.AtomGroup, build: md.BuildDir,
                        warnings: list[str]) -> tuple[np.ndarray, str]:
    """Human-OPRM1 residue ids for the protein, or tleap's 1..N if they cannot be transferred.

    `02.11.00` treats a failed transfer as fatal, because every residue number it reports would
    otherwise be wrong. A movie is exploratory and still worth watching with the wrong labels, so
    here it degrades to a warning -- but the JSON records WHICH numbering the PDB carries, so a
    residue picked out of the movie is never quietly mis-cited.
    """
    try:
        return md.construct_residue_map(u, build.receptor_pdb, notes=warnings), "human OPRM1 (P35372)"
    except (OSError, ValueError) as exc:
        warnings.append(f"could not transfer the construct numbering ({exc}); the movie PDB "
                        "carries the topology's own 1..N residue ids")
        return np.array([int(r.resid) for r in protein.residues], dtype=int), "topology (tleap 1..N)"


def write_topology(path: Path, subset: mda.AtomGroup, resids: np.ndarray,
                   groups: dict[str, mda.AtomGroup], slices: dict[str, list[int]],
                   numbering: str) -> dict[str, object]:
    """Write the movie's topology PDB: frame 0, with chains, HETATM flags and CONECT records.

    Three things a bare coordinate dump gets wrong, all of which a viewer then shows wrongly:

    * **Chains.** tleap concatenates every protein chain into one numbering with no TER, so a
      cartoon is drawn straight from the receptor's C-terminus into Gai's N-terminus as a bond
      across the cell. Each polypeptide is its own fragment (chains are not covalently joined),
      so the chains can be recovered from the bonding and written as separate chain IDs.
    * **ATOM vs HETATM.** The ligand, lipids, waters and ions are HETATM; left as ATOM, MolStar
      reads the macrocycle as part of the polymer and puts it in the cartoon.
    * **Residue numbering.** The protein is renumbered to the construct here, so a residue clicked
      in the movie is the residue the manuscript names.

    Returns the chain map, which goes into the JSON: both renderers need to know which chain is
    the RECEPTOR (drawn as the subject) and which are context, and the answer is the same question
    the pocket analysis asks -- the chain carrying the key residues -- so it is settled once here
    rather than re-guessed in matplotlib and again in JavaScript.
    """
    u = subset.universe
    for attr in ("chainIDs", "record_types", "elements"):
        if not hasattr(u.atoms, attr):
            u.add_TopologyAttr(attr)

    chains: dict[str, object] = {}
    fragments = list(groups["protein"].fragments)
    for i, frag in enumerate(fragments):
        frag.atoms.chainIDs = PROTEIN_CHAINS[i % len(PROTEIN_CHAINS)]
    for name, letter in (("ligand", CHAIN_LIGAND), ("lipid", CHAIN_LIPID),
                         ("water", CHAIN_WATER), ("ion", CHAIN_ION)):
        if len(groups[name]):
            groups[name].atoms.chainIDs = letter
            chains[name] = letter
    subset.atoms.record_types = "ATOM"
    for name in ("ligand", "lipid", "phosphate", "water", "ion"):
        start, stop = slices[name]
        if stop > start:
            subset[start:stop].record_types = "HETATM"
    if numbering.startswith("human"):
        groups["protein"].residues.resids = resids

    # Which protein chain is the receptor: the one carrying the pocket and functional residues.
    # Position in the file is not a safe proxy -- tleap's chain order follows the packing input,
    # not the biology -- and in system A the Gi heterotrimer supplies three more chains that must
    # be drawn as context rather than as the subject.
    protein_letters = [PROTEIN_CHAINS[i % len(PROTEIN_CHAINS)] for i in range(len(fragments))]
    receptor = sorted({letter for letter, frag in zip(protein_letters, fragments, strict=True)
                       if {int(r.resid) for r in frag.atoms.residues} & set(md.KEY_RESIDUES)})
    if not receptor and protein_letters:
        receptor = [protein_letters[0]]
    chains["receptor"] = receptor
    chains["other_protein"] = [c for c in protein_letters if c not in receptor]

    # bonds="all" rather than "conect": MDAnalysis writes only NON-guessed bonds under "conect",
    # which would silently drop every bond on the --topology smoke-test path where they had to be
    # guessed. Both sources are wanted here -- the point is that the viewer draws the macrocycle.
    with mda.coordinates.PDB.PDBWriter(str(path), bonds="all", reindex=True) as W:
        W.write(subset)
    return chains


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", action="append", default=[],
                    help="build directory (repeatable); default is every build found")
    ap.add_argument("--replica", help="export only this replica, e.g. prod_r1")
    ap.add_argument("--all", action="store_true", help="every replica of every build")
    ap.add_argument("--index", type=int, help="export the Nth (build, replica) pair, 1-based "
                                              "-- this is what the SLURM array task passes")
    ap.add_argument("--list", action="store_true", help="print the numbered pairs and exit")
    ap.add_argument("--frames", type=int, default=MOVIE_FRAMES,
                    help=f"frames in the movie (default {MOVIE_FRAMES}); the stride follows")
    ap.add_argument("--scan-frames", type=int, default=SCAN_FRAMES,
                    help="frames sampled to choose the lipid/water selection")
    ap.add_argument("--lipid-cut", type=float, default=LIPID_CUT,
                    help=f"A from the protein for an annular lipid (default {LIPID_CUT})")
    ap.add_argument("--water-cut", type=float, default=WATER_CUT,
                    help=f"A from the ligand/pocket for a kept water (default {WATER_CUT})")
    ap.add_argument("--waters", type=int, default=N_WATERS,
                    help=f"most-persistent pocket waters to keep (default {N_WATERS}); 0 for none")
    ap.add_argument("--no-lipids", action="store_true",
                    help="drop the annular lipids AND the phosphate reference plane: smallest "
                         "file, cleanest view of the receptor alone (waters are kept)")
    ap.add_argument("--hydrogens", action="store_true",
                    help="keep protein/ligand hydrogens (doubles the size; invisible in cartoon)")
    ap.add_argument("--no-align", action="store_true",
                    help="skip superposition -- keeps whole-box drift, tilt and PBC artefacts, "
                         "which is what you want when THOSE are the question")
    ap.add_argument("--topology", type=Path,
                    help="topology to use instead of the build's system.prmtop")
    ap.add_argument("--out", type=Path, default=OUT_ROOT)
    ap.add_argument("--force", action="store_true", help="re-export even if outputs exist")
    ap.add_argument("--all-builds", action="store_true",
                    help="include superseded builds, not just the newest per system")
    args = ap.parse_args()

    if args.build:
        builds = []
        for path in (Path(b).resolve() for b in args.build):
            found = md.discover_builds(path.parent, newest_only=False)
            match = [b for b in found if b.path == path]
            if not match:
                raise SystemExit(f"ERROR: {path} is not a build directory "
                                 "(expected <ligand>_<ASP|ASH>_<timestamp>/system.prmtop).")
            builds += match
    else:
        if not BUILD_ROOT.is_dir():
            raise SystemExit(f"ERROR: {BUILD_ROOT} does not exist. Build the systems first "
                             "(src/02.10.00_slurm_bundle/README.md).")
        builds = md.discover_builds(BUILD_ROOT, newest_only=not args.all_builds)

    todo = jobs(builds)
    if args.replica:
        todo = [j for j in todo if j.replica == args.replica]
    if not todo:
        raise SystemExit("ERROR: no production trajectories (prod_r*.dcd) found.")

    if args.list:
        for row in listing(todo, args.frames):
            print(row)
        for note in inventory_notes(builds):   # stderr: the listing stays one line per job
            print(note, file=sys.stderr)
        return 0

    if args.index is not None:
        if not 1 <= args.index <= len(todo):
            raise SystemExit(f"ERROR: --index {args.index} outside 1-{len(todo)}.")
        todo = [todo[args.index - 1]]
    elif not (args.all or args.build or args.replica):
        raise SystemExit("ERROR: nothing selected. Pass --all, --build, --index or --list.")

    failures = 0
    for job in todo:
        out_dir = args.out / job.build.name
        target = out_dir / f"{job.replica}_movie.json"
        if target.exists() and not args.force:
            print(f"skip {job.build.name}/{job.replica}: {target.name} exists (--force to redo)")
            continue
        print(f"--- {job.build.name} / {job.replica}: {job.dcd}", flush=True)
        try:
            summary = export_replica(
                job, out_dir, frames=args.frames, scan_frames=args.scan_frames,
                lipid_cut=0.0 if args.no_lipids else args.lipid_cut,
                water_cut=args.water_cut, n_waters=args.waters,
                hydrogens=args.hydrogens, membrane=not args.no_lipids,
                align=not args.no_align, topology=args.topology)
        except Exception as exc:  # noqa: BLE001
            # Deliberately broad, as in 02.11.00: a replica killed at its wall-time leaves a
            # TRUNCATED DCD and the readers raise whatever their parser happens to raise. In a job
            # array one unusable replica must not take the others with it.
            print(f"FAILED {job.build.name}/{job.replica}: {type(exc).__name__}: {exc}",
                  file=sys.stderr)
            failures += 1
            continue
        mv, at = summary["movie"], summary["atoms"]
        assert isinstance(mv, dict) and isinstance(at, dict)
        size_mb = (out_dir / f"{job.replica}_movie.xtc").stat().st_size / 1e6
        print(f"    {mv['n_frames']} frames every {mv['frame_spacing_ns']:.2f} ns, "
              f"{at['total']} atoms -> {size_mb:.1f} MB; {summary['wall_seconds']:.0f} s")
        for w in summary["warnings"]:  # type: ignore[union-attr]
            print(f"    WARNING: {w}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

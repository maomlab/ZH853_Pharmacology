#!/bin/bash
# Build the membrane system: PACKMOL-Memgen bilayer + tleap assembly -> Amber prmtop/rst7.
# Produces System A (receptor + ZH853 + Gi in POPC:chol 9:1). Run on a CPU node.
#
# IMPORTANT: do NOT `module load amber` here -- a system AmberTools conflicts with the conda
# env's AMBERHOME and breaks tool discovery (packmol-memgen: "reduce not available").
#   module unload amber cuda cudnn 2>/dev/null || true
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate zh853mor-prep   # AmberTools / PACKMOL-Memgen / obabel / reduce live here

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"

# --- build in a PRISTINE directory --------------------------------------------------------
# packmol-memgen writes to CWD *and reuses whatever it finds there*: the component PDBs
# (POPC/CHL1/WAT/Na+/Cl-, cf. its "Using WAT.pdb in the folder" message) and the preprocessed
# protein files (receptor_Trim_H.pdb, *.grid.pdb, receptorin_EMBED*.pdb). A stale component PDB
# whose atom order no longer matches the "atoms 1 20 / atoms 88 131" head/tail constraints in the
# generated packmol.inp makes PACKMOL fail with
#   "unable to put the molecules in the desired regions even without considering distance
#    tolerances" / "Maximum violation of the restraints: ~26" (= one lipid length, 100% of copies),
# and no change of packmol-memgen flags will fix it. So: always start empty.
# Generated artefacts live under intermediate/, never in src/ (SPECIFICATION D-16).
D250="${D250:-ASP}"   # ASP = charged D2.50 (default); ASH = protonated variant (D-11)
case "$D250" in ASP|ASH) ;; *) echo "ERROR: D250 must be ASP or ASH, got '$D250'."; exit 1 ;; esac
BUILD="${BUILD_DIR:-$REPO/intermediate/02.10.00_build/${D250}_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$BUILD"
[ -z "$(ls -A "$BUILD")" ] || { echo "ERROR: build dir $BUILD is not empty."; exit 1; }
echo "Building in $BUILD (D2.50 = $D250)"

# --- stage inputs -------------------------------------------------------------------------
# Receptor: use the membrane-ORIENTED file (normal along z) so --preoriented is valid.
# Produced by `make prep-orient` (src/02.05.00_orient_receptor.py); already carries the ACE/NME
# caps and the named His tautomers from 02.03.00, so no downstream default can override them.
# The D2.50 variant is a pure residue rename -- same geometry, so it is applied here rather than
# by duplicating the whole prep/orient chain (D-11).
if [ "$D250" = "ASH" ]; then
  awk '{ if ((substr($0,1,4)=="ATOM") && substr($0,18,3)=="ASP" && (substr($0,23,4)+0)==116)
           print substr($0,1,17) "ASH" substr($0,21); else print }' \
    "$REPO/intermediate/02.05.00_oriented/receptorR_oriented.pdb" > "$BUILD/receptor.pdb"
  grep -qc " ASH R 116" "$BUILD/receptor.pdb" || { echo "ERROR: ASP116->ASH rename found no atoms."; exit 1; }
else
  cp "$REPO/intermediate/02.05.00_oriented/receptorR_oriented.pdb" "$BUILD/receptor.pdb"
fi
cp "$HERE/tleap.in" "$HERE/make_tleap.py" "$HERE/fix_caps.py" "$HERE/check_placement.py" \
   "$HERE/check_piercing.py" "$BUILD/"
# Stage the run scripts too. The sbatch files invoke `python 02_equilibrate.py` by bare name and
# SLURM starts a job in the *submission* directory, so the build dir must be self-contained:
# everything for steps 3-5 is submitted from here, beside the system.prmtop it reads.
cp "$HERE/02_equilibrate.py" "$HERE/03_production.py" "$HERE/04_analyze.py" \
   "$HERE/submit_equilibrate.sbatch" "$HERE/submit_production.sbatch" \
   "$HERE/submit.sh" "$HERE/check_gpu_env.sh" "$HERE/cluster.env.example" "$BUILD/"
# cluster.env is the single source of the SLURM account/partition/GPU/wall-time. Stage it if it
# exists so the build directory is self-contained; submit.sh also falls back to the bundle copy,
# so a build made before cluster.env was filled in still works once it is.
if [ -f "$HERE/cluster.env" ]; then
  cp "$HERE/cluster.env" "$BUILD/"
else
  echo "NOTE: $HERE/cluster.env does not exist yet -- create it before step 3:"
  echo "        cp $HERE/cluster.env.example $HERE/cluster.env && \$EDITOR $HERE/cluster.env"
fi

# The finalised receptor does NOT arrive with `git pull`: intermediate/ is gitignored, and 02.03.00
# needs openmm/pdbfixer, which the cluster's zh853mor-prep env deliberately does not carry. So a
# stale copy from an earlier run survives a pull and silently rebuilds the system with charged
# termini and default HIE tautomers. Refuse to build from one.
missing=""
grep -q "^ATOM.* ACE R" "$BUILD/receptor.pdb" || missing="$missing ACE-cap"
grep -q "^ATOM.* NME R" "$BUILD/receptor.pdb" || missing="$missing NME-cap"
grep -qE "^ATOM.* (HID|HIE|HIP) R" "$BUILD/receptor.pdb" || missing="$missing His-tautomer-names"
if [ -n "$missing" ]; then
  echo "ERROR: the staged receptor is stale -- missing:$missing"
  echo "  $REPO/intermediate/02.05.00_oriented/receptorR_oriented.pdb predates D-15 (ACE/NME caps"
  echo "  and named His tautomers). It is gitignored, so a git pull does not update it."
  echo "  Regenerate it in the LOCAL analysis env (needs openmm + pdbfixer):"
  echo "      make prep-receptor prep-orient"
  echo "  then copy intermediate/02.05.00_oriented/receptorR_oriented.pdb to this machine."
  exit 1
fi
# Ligand parameters come from ligand_resp/run_resp.sh (run that first):
cp "$HERE"/ligand_resp/ZH853.mol2 "$HERE"/ligand_resp/ZH853.frcmod "$BUILD/" 2>/dev/null || \
  echo "WARNING: ZH853.mol2/.frcmod not found -- run ligand_resp/run_resp.sh first."

cd "$BUILD"

# --- pack the solvated bilayer around the oriented receptor ---
# POPC:CHL 9:1, 0.15 M NaCl, >=15 A water pad, extra z for the Gi domain (SPECIFICATION D-4/D-10).
# No --parametrize: packmol-memgen only builds/solvates the box; tleap.in does the FF assignment.
# NB: do NOT add --verbose -- packmol-memgen 2025.1.29 crashes on it ("NameError: streamer");
# for diagnostics inspect the generated packmol.inp / packmol-memgen.json instead.
packmol-memgen \
  --pdb receptor.pdb \
  --lipids POPC:CHL1 --ratio 9:1 \
  --salt --salt_c Na+ --saltcon 0.15 \
  --dist 15 --dist_wat 17.5 \
  --preoriented

# packmol-memgen writes bilayer_<input>.pdb; normalize the name for tleap.in.
# It exits 0 even when PACKMOL fails ("CRITICAL: No output file generated by PACKMOL"), so check.
if [ ! -s bilayer_receptor.pdb ]; then
  echo "ERROR: packmol-memgen produced no bilayer_receptor.pdb."
  echo "  Look for 'Maximum violation of the restraints' in packmol.log: a restraint-only failure"
  echo "  means the head/tail atom-index constraints in packmol.inp do not match the component PDBs,"
  echo "  not that the box is overpacked. Compare the staged POPC/CHL1/WAT PDBs against the"
  echo "  packmol_memgen library copies before touching --dist/--ratio."
  ls -1 bilayer_*.pdb 2>/dev/null
  exit 1
fi
mv -f bilayer_receptor.pdb bilayer_system.pdb

# memgen's preprocessing mangles the neutral caps: reduce leaves stray H on the capping amide
# (an HN2 survived the H-strip on NME), and the ff19SB libraries name the NME methyl carbon `C`
# while 02.03.00 writes the older pdbfixer name `CH3` -- loadpdb matches atoms by NAME, so
# tleap dies with "Atom .R<NME ...>.A<CH3> does not have a type". Normalise both caps to the
# installed library's template (renames and stray-H removal only; geometry is untouched).
python fix_caps.py bilayer_system.pdb receptor.pdb

# --- verify the OPM membrane registration survived the build (SPECIFICATION D-14) ------------
# Measured 2026-07-26: --preoriented is honoured, translation (+0.02, -0.18, +0.00) A, so the OPM
# frame from 02.05.00 survives. This is a regression guard for the case where it does not -- memgen
# re-centres the solute on its z bounding box when it orients the protein itself, and that centre is
# ~5 A below the OPM midplane for this receptor, which would embed it too high in the bilayer with
# nothing downstream complaining. check_placement.py measures the offset against the lipid phosphate
# planes and fails the build rather than passing a mis-embedded system to tleap. Set
# SKIP_PLACEMENT_CHECK=1 only to inspect a known-bad build.
if [ "${SKIP_PLACEMENT_CHECK:-0}" != "1" ]; then
  python check_placement.py bilayer_system.pdb receptor.pdb
fi

# --- lipid ring piercing --------------------------------------------------------------------
# packmol-memgen's own piercing finder is unreliable ("Lipid piercing finder failed" on this
# system). A tail threaded through an aromatic or sterol ring is topologically trapped: it cannot
# escape during minimisation or MD, so it silently corrupts the whole trajectory. Report only --
# a piercing needs a human decision (re-pack with a different seed vs edit the offender).
python check_piercing.py bilayer_system.pdb || true

# --- assemble with tleap (ff19SB + Lipid21 + OPC + GAFF2/RESP ligand) ------------------------
# make_tleap.py fills in the disulfide bond (loadpdb renumbers residues sequentially, so the
# OPRM1 numbering is gone) and the periodic box (must be PACKMOL's cell, not `setBox vdw`).
python make_tleap.py bilayer_system.pdb
tleap -f tleap_run.in

echo "Built system.prmtop / system.rst7 in $BUILD."
echo "Duplicate with D2.50 (Asp116) protonated for the parallel run."
cat <<EOF

Next (steps 3-5 run from the build directory, not from src/):
    cd $BUILD
    ./submit.sh check     # optional GPU pre-flight
    ./submit.sh all       # equilibration, then production chained with --dependency=afterok
Cluster account/partition/GPU/wall-time all come from cluster.env; edit that one file, not the
.sbatch scripts. \`./submit.sh all -n\` prints the sbatch commands without submitting.
EOF

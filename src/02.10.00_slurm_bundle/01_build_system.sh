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

# --- stage inputs into the current working directory (packmol-memgen writes to CWD) ---
# Receptor: use the membrane-ORIENTED file (normal along z) so --preoriented is valid.
# Produced by `make prep-orient` (src/02.05.00_orient_receptor.py). PPM/OPM is preferred for
# production; swap in the PPM-oriented PDB here if you have it.
cp "$REPO/intermediate/02.05.00_oriented/receptorR_oriented.pdb" ./receptor.pdb
# Ligand parameters come from ligand_resp/run_resp.sh (run that first):
cp ligand_resp/ZH853.mol2 ligand_resp/ZH853.frcmod . 2>/dev/null || \
  echo "WARNING: ZH853.mol2/.frcmod not found -- run ligand_resp/run_resp.sh first."

# --- pack the solvated bilayer around the oriented receptor ---
# POPC:CHL 9:1, 0.15 M NaCl, >=15 A water pad, extra z for the Gi domain (SPECIFICATION D-4/D-10).
# No --parametrize: packmol-memgen only builds/solvates the box; tleap.in does the FF assignment.
packmol-memgen \
  --pdb receptor.pdb \
  --lipids POPC:CHL1 --ratio 9:1 \
  --salt --salt_c Na+ --saltcon 0.15 \
  --dist 15 --dist_wat 17.5 \
  --preoriented

# packmol-memgen writes bilayer_<input>.pdb; normalize the name for tleap.in.
mv -f bilayer_receptor.pdb bilayer_system.pdb 2>/dev/null || \
  { echo "Check the packmol-memgen output name and update tleap.in accordingly:"; ls -1 bilayer_*.pdb; }

# --- assemble with tleap (ff19SB + Lipid21 + OPC + GAFF2/RESP ligand). See tleap.in. ---
tleap -f tleap.in

echo "Built system.prmtop / system.rst7. Duplicate with D2.50 (Asp116) protonated for the parallel run."

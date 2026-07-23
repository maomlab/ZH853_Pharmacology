#!/bin/bash
# Build the membrane system: PACKMOL-Memgen bilayer + tleap assembly -> Amber prmtop/rst7.
# Produces System A (receptor + ZH853 + Gi in POPC:chol 9:1). Run on a CPU node.
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate zh853mor-sim

REC=receptorR_fixed_heavy.pdb          # from intermediate/02.03.00_receptor/
LIG_MOL2=ZH853.mol2                     # from ligand_resp/ (GAFF2/RESP charges)
LIG_FRC=ZH853.frcmod

# 1) Pre-orient to the membrane normal (PPM/OPM). TODO: fetch PPM output or use the cholesterol
#    frame computed locally (product/02.01.00_prep_assessment_*: normal + span).

# 2) Pack the bilayer around the (pre-oriented) receptor. POPC:CHL 9:1, 0.15 M NaCl, OPC water,
#    >=15 A water pad, extra intracellular Z for the Gi domain (SPECIFICATION D-4/D-10).
packmol-memgen \
  --pdb "${REC}" \
  --lipids POPC:CHL1 --ratio 9:1 \
  --salt --salt_c Na+ --saltcon 0.15 \
  --dist 15 --dist_wat 17.5 \
  --preoriented \
  --parametrize   # emits Lipid21/ff19SB leap; we override ligand+water below

# 3) Assemble with tleap (ff19SB + Lipid21 + OPC + GAFF2 ligand). See tleap.in.
tleap -f tleap.in

echo "Built system.prmtop / system.rst7. Duplicate with D2.50 (Asp116) protonated for the parallel run."

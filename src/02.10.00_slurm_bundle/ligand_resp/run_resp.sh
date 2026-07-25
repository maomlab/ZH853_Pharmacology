#!/bin/bash
# ZH853 ligand parameterization. Route A (quick) is the default; Route B (RESP) for FEP.
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate zh853mor-prep   # antechamber/parmchk2/RESP live in the prep env

SDF=../../../intermediate/02.04.00_ligand/ZH853_prepared.sdf   # +1, explicit H, bond orders

# ---- Route A: GAFF2 + AM1-BCC (minutes; good for equilibration/MD) ----
antechamber -i "${SDF}" -fi sdf -o ZH853.mol2 -fo mol2 -c bcc -nc 1 -at gaff2 -rn LIG
parmchk2 -i ZH853.mol2 -f mol2 -o ZH853.frcmod -s gaff2
echo "Route A done: ZH853.mol2 + ZH853.frcmod"

# ---- Route B: multi-conformer RESP (HF/6-31G*), for FEP charge accuracy (SPECIFICATION D-3) ----
# 1) generate low-energy conformers (macrocycle-aware) from the SDF
# 2) QM ESP per conformer (Psi4 or Gaussian):  # TODO(OQ-3): which QM engine is available?
#      psi4 esp.dat  (HF/6-31G* single points on ESP grid)
# 3) two-stage RESP fit over conformers:
#      antechamber ... -c resp   /   resp -O -i resp1.in ...
# Uncomment and wire once the QM engine + queue are confirmed.

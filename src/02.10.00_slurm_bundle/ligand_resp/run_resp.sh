#!/bin/bash
# Ligand parameterization for one of the four ZH cyclic peptides.
# Route A (GAFF2/AM1-BCC, minutes) is the default; Route B (RESP) is for FEP charge accuracy.
#
#   ./ligand_resp/run_resp.sh ZH853        # one ligand
#   for L in ZH853 ZH850 ZH831 ZH809; do ./ligand_resp/run_resp.sh $L; done
#
# Output goes to ligand_resp/<LIGAND>/<LIGAND>.{mol2,frcmod}, which 01_build_system.sh stages
# into the build directory. Per-ligand subdirectories because antechamber writes a pile of
# fixed-name scratch files (ANTECHAMBER_*.AC, sqm.in, sqm.out) that would collide otherwise.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BUNDLE="$(cd "$HERE/.." && pwd)"
REPO="$(cd "$BUNDLE/../.." && pwd)"

LIGAND="${1:-}"
if [ -z "$LIGAND" ]; then
  echo "usage: $0 <ligand>   (one of: $(python "$BUNDLE/ligands.py" --list | tr -d '\n'))" >&2
  exit 1
fi
if [ "$LIGAND" = "apo" ]; then
  echo "ERROR: 'apo' has no ligand to parameterize." >&2
  exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ZH_PREP_ENV:-zh853mor-prep}"   # antechamber/parmchk2/RESP live in the prep env

SDF="$REPO/$(python "$BUNDLE/ligands.py" --field sdf --ligand "$LIGAND")"
CHG="$(python "$BUNDLE/ligands.py" --field net_charge --ligand "$LIGAND")"
if [ ! -f "$SDF" ]; then
  echo "ERROR: prepared SDF not found: $SDF" >&2
  echo "  Generate it in the LOCAL analysis env first (needs RDKit):  make prep-ligand" >&2
  exit 1
fi

OUT="$HERE/$LIGAND"
mkdir -p "$OUT"
cd "$OUT"

# ---- Route A: GAFF2 + AM1-BCC (minutes; good for equilibration/MD) ----
# -rn LIG: the unit name tleap matches the packed residue against. fix_ligand.py renames the
# packed ligand (deposited name L01) to LIG to meet it. Do NOT change one without the other.
echo "Parameterizing $LIGAND (net charge $CHG) from $(basename "$SDF")"
antechamber -i "$SDF" -fi sdf -o "$LIGAND.mol2" -fo mol2 -c bcc -nc "$CHG" -at gaff2 -rn LIG
parmchk2 -i "$LIGAND.mol2" -f mol2 -o "$LIGAND.frcmod" -s gaff2

# parmchk2 emits ATTN lines for parameters it had to guess; those are the ones worth reading.
if grep -q "ATTN" "$LIGAND.frcmod"; then
  echo "NOTE: $LIGAND.frcmod contains ATTN (guessed) parameters:"
  grep -n "ATTN" "$LIGAND.frcmod" | head -10
fi
echo "Route A done: $OUT/$LIGAND.mol2 + $LIGAND.frcmod"

# ---- Route B: multi-conformer RESP (HF/6-31G*), for FEP charge accuracy (SPECIFICATION D-3) ----
# 1) generate low-energy conformers (macrocycle-aware) from the SDF
# 2) QM ESP per conformer (Psi4 or Gaussian):  # TODO(OQ-3): which QM engine is available?
#      psi4 esp.dat  (HF/6-31G* single points on ESP grid)
# 3) two-stage RESP fit over conformers:
#      antechamber ... -c resp   /   resp -O -i resp1.in ...
# Uncomment and wire once the QM engine + queue are confirmed.

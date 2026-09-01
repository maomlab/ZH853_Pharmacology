#!/bin/bash
# GAFF2 + AM1-BCC force-field parameters for the ZH ligands (Phase 2).
#
#   bash src/02.08.00_ligand_parameterize.sh                 # all ligands in the registry
#   bash src/02.08.00_ligand_parameterize.sh ZH853 ZH850     # a subset
#   make prep-ligand-parameterize
#
# Writes intermediate/02.08.00_ligand_params/<LIGAND>/<LIGAND>.{mol2,frcmod}, which
# src/02.10.00_slurm_bundle/01_build_system.sh stages into each build directory. Per-ligand
# subdirectories because antechamber writes fixed-name scratch (ANTECHAMBER_*.AC, sqm.in,
# sqm.out) that would otherwise collide between ligands.
#
# Was src/02.10.00_slurm_bundle/ligand_resp/run_resp.sh. Moved out of the bundle because it needs
# neither a GPU nor SLURM, and renamed because only the (still stubbed) Route B uses RESP -- what
# actually runs today is AM1-BCC. Outputs now land in intermediate/, not in src/: src/ is code.
#
# ENVIRONMENT: needs antechamber + parmchk2 from AmberTools, which are in zh853mor-prep, NOT in
# the local analysis env. This is the one Phase-2 prep step that does not run in zh853mor-local,
# which is also why it is not part of `make prep`.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
# Single registry of ligand names, net charges and input paths. It lives in the bundle because the
# bundle must stay self-contained on the cluster; querying it here avoids a second copy that could
# drift from it.
REGISTRY="$REPO/src/02.10.00_slurm_bundle/ligands.py"
OUTROOT="$REPO/intermediate/02.08.00_ligand_params"

for tool in antechamber parmchk2; do
  command -v "$tool" >/dev/null 2>&1 || {
    echo "ERROR: $tool not found on PATH." >&2
    echo "  It comes from AmberTools, which is in the zh853mor-prep env, not zh853mor-local:" >&2
    echo "      conda activate zh853mor-prep" >&2
    echo "  (create it with \`make env-cluster\`; see src/02.10.00_slurm_bundle/README.md)" >&2
    exit 1
  }
done

LIGANDS=("$@")
if [ ${#LIGANDS[@]} -eq 0 ]; then
  # every registry entry except the apo pseudo-system, which has no ligand to parameterize
  read -r -a LIGANDS <<< "$(python "$REGISTRY" --list | sed 's/\bapo\b//')"
fi

for LIGAND in "${LIGANDS[@]}"; do
  [ -n "$LIGAND" ] || continue
  if [ "$LIGAND" = "apo" ]; then
    echo "skipping 'apo': no ligand to parameterize"
    continue
  fi

  SDF="$REPO/$(python "$REGISTRY" --field sdf --ligand "$LIGAND")"
  CHG="$(python "$REGISTRY" --field net_charge --ligand "$LIGAND")"
  if [ ! -f "$SDF" ]; then
    echo "ERROR: prepared SDF not found: $SDF" >&2
    echo "  Generate it first:  make prep-ZH853-protonate   (ZH853)" >&2
    echo "                      make prep-analogs-pose      (ZH850/ZH831/ZH809)" >&2
    exit 1
  fi

  OUT="$OUTROOT/$LIGAND"
  mkdir -p "$OUT"
  echo "Parameterizing $LIGAND (net charge $CHG) from $(basename "$SDF") -> ${OUT#"$REPO"/}"

  # -rn LIG: the unit name tleap matches the packed residue against. fix_ligand.py renames the
  # packed ligand (deposited name L01) to LIG to meet it. Do NOT change one without the other.
  (
    cd "$OUT"
    antechamber -i "$SDF" -fi sdf -o "$LIGAND.mol2" -fo mol2 -c bcc -nc "$CHG" -at gaff2 -rn LIG
    parmchk2 -i "$LIGAND.mol2" -f mol2 -o "$LIGAND.frcmod" -s gaff2
  ) 2>&1 | tee "$OUT/antechamber.log"

  # Do not trust exit status alone: AmberTools tools have been known to exit 0 having written
  # nothing useful, and an empty mol2 would surface much later as an unparameterized residue in
  # tleap. Same reason 01_build_system.sh checks for packmol-memgen's output file.
  for f in "$LIGAND.mol2" "$LIGAND.frcmod"; do
    [ -s "$OUT/$f" ] || {
      echo "ERROR: $OUT/$f was not produced (or is empty)." >&2
      echo "  Read $OUT/antechamber.log; sqm convergence failures show up there." >&2
      exit 1
    }
  done

  # parmchk2 emits ATTN for parameters it had to guess; those are the ones worth reading, and the
  # ones to revisit if a ligand behaves oddly in MD.
  if grep -q "ATTN" "$OUT/$LIGAND.frcmod"; then
    echo "  NOTE: $LIGAND.frcmod contains guessed (ATTN) parameters:"
    grep -n "ATTN" "$OUT/$LIGAND.frcmod" | head -10 | sed 's/^/    /'
  fi
  echo "  done: $LIGAND.mol2 + $LIGAND.frcmod + antechamber.log"
done

cat <<EOM

Route B (multi-conformer RESP, HF/6-31G*) is not wired up: it needs a QM engine
(Psi4 or Gaussian) that has not been chosen -- the remaining TODO(OQ-3). AM1-BCC above is
appropriate for equilibration and production MD; RESP matters for FEP charge accuracy
(SPECIFICATION D-3).

Next: build the membrane systems (src/02.10.00_slurm_bundle, on the cluster):
    cd src/02.10.00_slurm_bundle && LIGAND=ZH853 D250=ASP ./01_build_system.sh
EOM

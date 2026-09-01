#!/bin/bash
#SBATCH --job-name=zh853_build
#SBATCH --output=build_%A_%a.out
#
# CPU stage: build the membrane systems, ONE (LIGAND, D2.50) COMBINATION PER ARRAY TASK, so the
# ten builds run concurrently instead of ~2.5 h of PACKMOL-Memgen serially on a login node.
#
#     ./submit.sh build
#
# Parallel builds are safe by construction: each writes to its own timestamped directory, which is
# exactly why 01_build_system.sh insists on a pristine one -- PACKMOL-Memgen writes to the CWD and
# silently reuses whatever it finds there. Directory names carry LIGAND and D250, so tasks cannot
# collide even when they start in the same second.

set -euo pipefail

# shellcheck source=cluster_env.sh
source "${SLURM_SUBMIT_DIR:-$PWD}/cluster_env.sh" || exit 1

BUNDLE="${SLURM_SUBMIT_DIR:-$PWD}"

# The (ligand x D2.50) grid, in registry order, so an array index always means the same build.
COMBOS=()
for L in $(python "$BUNDLE/ligands.py" --list); do
  for D in ASP ASH; do COMBOS+=("$L:$D"); done
done
IDX="${SLURM_ARRAY_TASK_ID:-1}"
COMBO="${COMBOS[$((IDX - 1))]:-}"
[ -n "$COMBO" ] || { echo "ERROR: no combination at array index $IDX (have ${#COMBOS[@]})."; exit 1; }
LIGAND="${COMBO%%:*}"
D250="${COMBO##*:}"

echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  task ${IDX}/${#COMBOS[@]}  ${LIGAND} / ${D250}"
cd "$BUNDLE"
LIGAND="$LIGAND" D250="$D250" ./01_build_system.sh

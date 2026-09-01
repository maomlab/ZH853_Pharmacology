#!/bin/bash
#SBATCH --job-name=zh853_preprod
#SBATCH --output=preprod_%j.out
#
# Step 3.5: UNRESTRAINED pre-production. Submit with the wrapper, from the BUILD DIRECTORY,
# after step 3:
#
#     ./submit.sh preprod       # or: ./submit.sh all   (chains eq -> preprod -> prod)
#
# Why this exists: the six-stage ramp in step 3 is ~2.25 ns, adapted from CHARMM-GUI -- which
# hands its schedule a PRE-EQUILIBRATED bilayer from its own builder. Ours comes straight from
# PACKMOL, where lipids are packed in an artificially ordered, extended conformation. The protein
# relaxes in a few hundred ps; area per lipid takes tens of ns. This leg is that time, and it is
# DISCARDED -- it is equilibration, not sampling. Production then resumes from preprod_final.xml.
#
# It also fixes a smaller problem: stage 6 of the ramp still holds the backbone at
# 10 kJ/mol/nm^2, so without this leg production would start from a state never sampled unrestrained.
#
# Physics is identical to production (03_production.py: NPT, membrane barostat, HMR/4 fs), so
# there is no separate script -- only the length, the seed and the discard.

set -euo pipefail

# cluster.env is REQUIRED: running on built-in defaults would use a different conda env and a
# different system basename and still exit 0. cluster_env.sh fails loudly instead.
# shellcheck source=cluster_env.sh
source "${SLURM_SUBMIT_DIR:-$PWD}/cluster_env.sh" || exit 1

module purge
if [ -n "${ZH_MODULES:-}" ]; then
  # shellcheck disable=SC2086 -- ZH_MODULES is a space-separated module list
  module load ${ZH_MODULES}
fi

source "${ZH_CONDA_SH:-$(conda info --base)/etc/profile.d/conda.sh}"
conda activate "${ZH_SIM_ENV:-zh853mor-sim}"

SYS="${ZH_SYS:-system}"
NS="${ZH_PREPROD_NS:-100}"

# Pre-production RESUMES rather than restarts. The membrane is the slow degree of freedom and one
# leg is often not enough -- but re-running this script used to start again from the restrained
# ramp, silently discarding the leg already computed. Legs are preprod, preprod2, preprod3 ...,
# each continuing from the previous one's final state, so "it needs longer" costs only the extra
# time. Re-run this stage as many times as the QC says.
LEG=1
PREV=""
while :; do
  NAME=$([ "$LEG" -eq 1 ] && echo preprod || echo "preprod${LEG}")
  [ -f "${NAME}_final.xml" ] || break
  PREV="$NAME"
  LEG=$((LEG + 1))
done
OUT=$([ "$LEG" -eq 1 ] && echo preprod || echo "preprod${LEG}")
if [ -n "$PREV" ]; then
  STATE="${PREV}_final.xml"
  echo "continuing pre-production from ${STATE} (leg ${LEG}); cumulative >= $(( (LEG - 1) * NS + NS )) ns"
else
  STATE="${SYS}_eq.xml"
  echo "starting pre-production from the restrained ramp (${STATE})"
fi
[ -f "$STATE" ] || { echo "ERROR: state file $STATE not found in $PWD."; exit 1; }

echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  unrestrained: ${NS} ns -> ${OUT}"
python -c "import openmm; print('openmm', openmm.version.version)"

# seed 0 marks the discard leg; production replicas use seeds 1..ZH_REPLICAS.
python 03_production.py \
  --prmtop "${SYS}.prmtop" \
  --state  "${STATE}" \
  --out    "${OUT}" \
  --ns "${NS}" --seed 0

# Re-assess on the equilibrated-for-real trajectory. With ZH_QC_GATE=1 (default) a FAIL exits
# non-zero, so a production job chained with --dependency=afterok never launches on a broken
# system -- the whole point of checking before spending replica-days of GPU time. A WARN exits 0.
echo
if [ "${ZH_QC_GATE:-1}" = "1" ]; then
  python check_equilibration.py --top "${SYS}.prmtop" --eq "${OUT}" --receptor receptor.pdb --out "${OUT}_qc"
else
  python check_equilibration.py --top "${SYS}.prmtop" --eq "${OUT}" --receptor receptor.pdb --out "${OUT}_qc" || true
fi

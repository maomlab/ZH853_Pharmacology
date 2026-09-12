#!/bin/bash
# Submit the trajectory reduction as a job array, one production replica per task.
#
#     ./submit_reduce.sh            # every replica of every (newest) build
#     ./submit_reduce.sh -n         # dry run: print the sbatch command, submit nothing
#     ./submit_reduce.sh -- --force # anything after -- is passed to 01_reduce_trajectory.py
#
# Resources come from cluster.env at the repository ROOT (the CPU settings, ZH_CPU_*), the same
# file the build and simulation stages use. The array size is taken from the script's own listing
# of (build, replica) pairs, so adding a replica or a system needs no edit here.
#
# Reduction is CPU- and I/O-bound: ~5 GB read per 500 ns replica, a few minutes of compute.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
die() { echo "ERROR: $*" >&2; exit 1; }

DRY=0
if [ "${1:-}" = "-n" ]; then DRY=1; shift; fi
[ "${1:-}" = "--" ] && shift
EXTRA=("$@")

# shellcheck source=../02.10.00_slurm_bundle/cluster_env.sh
source "$HERE/../02.10.00_slurm_bundle/cluster_env.sh" || exit 1
[ -n "${ZH_ACCOUNT:-}" ] || die "ZH_ACCOUNT is empty in ${ZH_CLUSTER_ENV}."
: "${ZH_CPU_TIME:=04:00:00}" ; : "${ZH_CPU_CPUS:=4}" ; : "${ZH_CPU_MEM:=16G}"
: "${ZH_CPU_PARTITION:=${ZH_GPU_PARTITION:-}}"
[ -n "$ZH_CPU_PARTITION" ] || die "set ZH_CPU_PARTITION (or ZH_GPU_PARTITION) in ${ZH_CLUSTER_ENV}."

# The listing is the source of truth for the array size AND for what each index means.
LISTING="$(python "$HERE/01_reduce_trajectory.py" --list)" \
  || die "could not list the production replicas (is intermediate/02.10.00_build populated?)."
N=$(printf '%s\n' "$LISTING" | grep -c .)
[ "$N" -gt 0 ] || die "no production trajectories found."
printf '%s\n' "$LISTING" | sed 's/^/  /'

SBATCH_ARGS=(
  --account="$ZH_ACCOUNT"
  --partition="$ZH_CPU_PARTITION"
  --time="$ZH_CPU_TIME"
  --cpus-per-task="$ZH_CPU_CPUS"
  --mem="$ZH_CPU_MEM"
  --array="1-${N}"
  --export="ALL,ZH_CLUSTER_ENV=$ZH_CLUSTER_ENV"
)
[ -n "${ZH_QOS:-}" ]        && SBATCH_ARGS+=(--qos="$ZH_QOS")
[ -n "${ZH_CONSTRAINT:-}" ] && SBATCH_ARGS+=(--constraint="$ZH_CONSTRAINT")
if [ -n "${ZH_EXTRA_SBATCH:-}" ]; then
  # shellcheck disable=SC2206 -- deliberate word-splitting: ZH_EXTRA_SBATCH holds whole flags
  SBATCH_ARGS+=($ZH_EXTRA_SBATCH)
fi

if [ "$DRY" -eq 1 ]; then
  echo "DRY RUN: sbatch ${SBATCH_ARGS[*]} $HERE/submit_reduce.sbatch ${EXTRA[*]-}"
  exit 0
fi
jid=$(cd "$HERE" && sbatch --parsable "${SBATCH_ARGS[@]}" submit_reduce.sbatch "${EXTRA[@]-}")
echo "submitted reduction: $jid   ($N replicas, one per array task)"
echo "  -> intermediate/02.11.00_analysis/<system>/<replica>.{npz,json}"
echo "then, on a machine with the LOCAL env:"
echo "    make sim-aggregate sim-figures"

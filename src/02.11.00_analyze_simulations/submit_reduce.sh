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
# Pass-through flags for 01_reduce_trajectory.py. Expanded below as ${EXTRA[@]+"${EXTRA[@]}"}:
# "${EXTRA[@]-}" looks equivalent but expands an EMPTY array to one empty string, which sbatch
# forwards to the job script and argparse then rejects with `unrecognized arguments:` -- every
# task of the array failing on an argument nobody passed.
EXTRA=("$@")

# shellcheck source=../02.10.00_slurm_bundle/cluster_env.sh
source "$HERE/../02.10.00_slurm_bundle/cluster_env.sh" || exit 1
[ -n "${ZH_ACCOUNT:-}" ] || die "ZH_ACCOUNT is empty in ${ZH_CLUSTER_ENV}."
: "${ZH_CPU_TIME:=04:00:00}" ; : "${ZH_CPU_CPUS:=4}" ; : "${ZH_CPU_MEM:=16G}"
: "${ZH_CPU_PARTITION:=${ZH_GPU_PARTITION:-}}"
[ -n "$ZH_CPU_PARTITION" ] || die "set ZH_CPU_PARTITION (or ZH_GPU_PARTITION) in ${ZH_CLUSTER_ENV}."

# The listing is the source of truth for the array size AND for what each index means. It also
# carries each replica's LENGTH, read from the DCD header and the state log: a replica that hit
# its wall-time or is still being written is indistinguishable from a finished one by file name
# alone, and reducing it silently analyses a partial run. Per-build remarks (chiefly missing
# replicas) come back on stderr, so they pass through to the terminal rather than into $LISTING.
LISTING="$(python "$HERE/01_reduce_trajectory.py" --list)" \
  || die "could not list the production replicas (is intermediate/02.10.00_build populated?)."
N=$(printf '%s\n' "$LISTING" | grep -c .)
[ "$N" -gt 0 ] || die "no production trajectories found."

FMT='  %-3s %-12s %-9s %8s %9s %9s  %s\n'
# shellcheck disable=SC2059 -- FMT is a fixed format string, not user input
printf "$FMT" "#" "system" "replica" "frames" "ns" "target" "status"
printf '%s\n' "$LISTING" | awk -F'\t' -v fmt="$FMT" \
  '{ printf fmt, $1, $2, $3, $4, $5, $6, $7 }'

read -r TOTAL_NS N_PARTIAL N_WRITING N_BAD <<SUMMARY
$(printf '%s\n' "$LISTING" | awk -F'\t' '
   { if ($5 != "?") total += $5
     if ($7 ~ /^partial/)          partial++
     else if ($7 == "writing")     writing++
     else if ($7 ~ /^unreadable/)  bad++ }
   END { printf "%.0f %d %d %d", total, partial + 0, writing + 0, bad + 0 }')
SUMMARY
echo "  ${N} replicas, ${TOTAL_NS} ns of trajectory in total"

if [ "$N_WRITING" -gt 0 ]; then
  echo "WARNING: $N_WRITING replica(s) were written within the last 15 minutes -- the production"
  echo "  job is probably still running. Reducing one now measures whatever it has reached so far."
fi
if [ "$N_PARTIAL" -gt 0 ]; then
  echo "WARNING: $N_PARTIAL replica(s) are short of the ZH_PROD_NS their build was configured for."
  echo "  A run that hit its wall-time can be extended; reducing it now analyses the partial run."
fi
if [ "$N_BAD" -gt 0 ]; then
  echo "WARNING: $N_BAD trajectory file(s) could not be read at all (truncated or still opening)."
fi
if [ $((N_WRITING + N_PARTIAL + N_BAD)) -gt 0 ]; then
  echo "  To reduce only what is finished, pass --build <dir> (repeatable) or --replica <name>"
  echo "  to 01_reduce_trajectory.py directly instead of submitting the whole array."
fi

# SLURM writes --output relative to the submission directory, which is this stage directory under
# src/ -- where generated files do not belong (D-16). Keep the logs with the data this stage
# produces instead. The directory must exist before sbatch runs: a job whose output file cannot be
# opened does not start.
LOGS="$(cd "$HERE/../.." && pwd)/intermediate/02.11.00_analyze_simulations/logs"
[ "$DRY" -eq 1 ] || mkdir -p "$LOGS" || die "cannot create the log directory $LOGS."

SBATCH_ARGS=(
  --account="$ZH_ACCOUNT"
  --output="$LOGS/reduce_%A_%a.out"
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
  echo "DRY RUN: sbatch ${SBATCH_ARGS[*]} $HERE/submit_reduce.sbatch" \
       "${EXTRA[@]+${EXTRA[*]}}"
  exit 0
fi
jid=$(cd "$HERE" && sbatch --parsable "${SBATCH_ARGS[@]}" submit_reduce.sbatch \
        ${EXTRA[@]+"${EXTRA[@]}"})
echo "submitted reduction: $jid   ($N replicas, one per array task)"
echo "  -> intermediate/02.11.00_analyze_simulations/<system>/<replica>.{npz,json}"
echo "  logs: intermediate/02.11.00_analyze_simulations/logs/reduce_${jid}_*.out"
echo "then, on a machine with the LOCAL env:"
echo "    make sim-aggregate sim-figures"

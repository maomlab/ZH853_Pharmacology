#!/bin/bash
# Submit the movie export as a job array, one production replica per task.
#
#     ./submit_export.sh              # every replica of every (newest) build
#     ./submit_export.sh -n           # dry run: print the sbatch command, submit nothing
#     ./submit_export.sh -- --force   # anything after -- is passed to 01_export_movie_trajectory.py
#     ./submit_export.sh -- --frames 600 --no-align
#
# Resources come from cluster.env at the repository ROOT (the CPU settings, ZH_CPU_*), the same
# file the build, simulation and reduction stages use. The array size is taken from the script's
# own listing of (build, replica) pairs, so adding a replica or a system needs no edit here.
#
# Export is CPU- and I/O-bound: the same ~5 GB read per 500 ns replica as the reduction, and a
# few minutes of imaging and superposition. The OUTPUT is what matters -- a few tens of MB per
# replica, which is what gets copied back to a laptop.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
die() { echo "ERROR: $*" >&2; exit 1; }

DRY=0
if [ "${1:-}" = "-n" ]; then DRY=1; shift; fi
[ "${1:-}" = "--" ] && shift
# Pass-through flags for 01_export_movie_trajectory.py. Expanded below as ${EXTRA[@]+"${EXTRA[@]}"}:
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

# The listing is the source of truth for the array size AND for what each index means. It carries
# each replica's LENGTH, read from the DCD header, and what the movie will be: how many frames
# survive the stride, and how much simulated time each of them spans. That last column is the one
# to read before submitting -- a 500 ns replica in 300 frames samples every 1.7 ns, and anything
# faster than that (a sidechain flip, a water exchange) is simply not in the result.
#
# --frames is forwarded from EXTRA if it was given, so the listing describes the movie that will
# actually be made rather than the default one.
FRAMES=300
for i in "${!EXTRA[@]}"; do
  [ "${EXTRA[$i]}" = "--frames" ] && FRAMES="${EXTRA[$((i + 1))]:-$FRAMES}"
done
LISTING="$(python "$HERE/01_export_movie_trajectory.py" --list --frames "$FRAMES")" \
  || die "could not list the production replicas (is intermediate/02.10.00_build populated?)."
N=$(printf '%s\n' "$LISTING" | grep -c .)
[ "$N" -gt 0 ] || die "no production trajectories found."

FMT='  %-3s %-12s %-9s %8s %9s %8s %10s %7s  %s\n'
# shellcheck disable=SC2059 -- FMT is a fixed format string, not user input
printf "$FMT" "#" "system" "replica" "frames" "ns" "movie" "ns/frame" "stride" "status"
printf '%s\n' "$LISTING" | awk -F'\t' -v fmt="$FMT" \
  '{ printf fmt, $1, $2, $3, $4, $5, $6, $7, $8, $9 }'

read -r TOTAL_NS N_PARTIAL N_WRITING N_BAD <<SUMMARY
$(printf '%s\n' "$LISTING" | awk -F'\t' '
   { if ($5 != "?") total += $5
     if ($9 ~ /^partial/)          partial++
     else if ($9 == "writing")     writing++
     else if ($9 ~ /^unreadable/)  bad++ }
   END { printf "%.0f %d %d %d", total, partial + 0, writing + 0, bad + 0 }')
SUMMARY
echo "  ${N} replicas, ${TOTAL_NS} ns of trajectory in total"

if [ "$N_WRITING" -gt 0 ]; then
  echo "WARNING: $N_WRITING replica(s) were written within the last 15 minutes -- the production"
  echo "  job is probably still running. A movie made now stops wherever the run has reached."
fi
if [ "$N_PARTIAL" -gt 0 ]; then
  echo "WARNING: $N_PARTIAL replica(s) are short of the ZH_PROD_NS their build was configured for."
  echo "  The movie's clock will end early; the caption reports the time it actually reaches."
fi
if [ "$N_BAD" -gt 0 ]; then
  echo "WARNING: $N_BAD trajectory file(s) could not be read at all (truncated or still opening)."
fi
if [ $((N_WRITING + N_PARTIAL + N_BAD)) -gt 0 ]; then
  echo "  To export only what is finished, pass --build <dir> (repeatable) or --replica <name>"
  echo "  to 01_export_movie_trajectory.py directly instead of submitting the whole array."
fi

# SLURM writes --output relative to the submission directory, which is this stage directory under
# src/ -- where generated files do not belong (D-16). Keep the logs with the data this stage
# produces instead. The directory must exist before sbatch runs: a job whose output file cannot be
# opened does not start.
LOGS="$(cd "$HERE/../.." && pwd)/intermediate/02.12.00_render_trajectory_movies/logs"
[ "$DRY" -eq 1 ] || mkdir -p "$LOGS" || die "cannot create the log directory $LOGS."

SBATCH_ARGS=(
  --account="$ZH_ACCOUNT"
  --output="$LOGS/export_%A_%a.out"
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
  echo "DRY RUN: sbatch ${SBATCH_ARGS[*]} $HERE/submit_export.sbatch" \
       "${EXTRA[@]+${EXTRA[*]}}"
  exit 0
fi
jid=$(cd "$HERE" && sbatch --parsable "${SBATCH_ARGS[@]}" submit_export.sbatch \
        ${EXTRA[@]+"${EXTRA[@]}"})
echo "submitted movie export: $jid   ($N replicas, one per array task)"
echo "  -> intermediate/02.12.00_render_trajectory_movies/<system>/<replica>_movie.{pdb,xtc,json}"
echo "  logs: intermediate/02.12.00_render_trajectory_movies/logs/export_${jid}_*.out"
echo "then copy that directory to the machine with the LOCAL env and:"
echo "    make movies"

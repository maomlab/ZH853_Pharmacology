#!/bin/bash
# Locate and load cluster.env.  *** SOURCE this file; do not execute it. ***
#
#     source "${SLURM_SUBMIT_DIR:-$PWD}/cluster_env.sh" || exit 1
#
# Every script that reads a ZH_* setting REQUIRES cluster.env. Silently falling back to built-in
# defaults is worse than stopping: the job would activate a different conda env, or run under a
# different system basename, or for a different number of ns -- and still exit 0 looking like it
# worked. A wrong answer that reports success is the one failure mode worth spending lines to
# prevent, which is why this is an error and not a warning.
#
# The exception is 01_build_system.sh: it is a CPU build step that needs no SLURM settings, so it
# sets ZH_ENV_OPTIONAL=1 and continues (with a note) when the file is absent. That is also why the
# file can be created after a build without invalidating it.
#
# Lookup order, first hit wins:
#   $ZH_CLUSTER_ENV     -- exported by submit.sh, so a job uses exactly what the submission used
#   $PWD                -- the build directory (SLURM starts a job in the submission directory)
#   $SLURM_SUBMIT_DIR   -- explicit fallback if something changed the CWD
#   this script's own directory
# Note that job scripts cannot find this file relative to themselves: sbatch copies the job script
# to a spool directory, so $0 there is not the build directory. $PWD is.

_zh_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)" || _zh_dir="$PWD"
_zh_found=""
for _zh_cand in "${ZH_CLUSTER_ENV:-}" "$PWD/cluster.env" \
                "${SLURM_SUBMIT_DIR:+$SLURM_SUBMIT_DIR/cluster.env}" "$_zh_dir/cluster.env"; do
  [ -n "$_zh_cand" ] || continue
  if [ -f "$_zh_cand" ]; then
    _zh_found="$(cd "$(dirname "$_zh_cand")" && pwd)/$(basename "$_zh_cand")"
    break
  fi
done

if [ -n "$_zh_found" ]; then
  # shellcheck disable=SC1090
  . "$_zh_found"
  ZH_CLUSTER_ENV="$_zh_found"
  export ZH_CLUSTER_ENV
  echo "cluster settings: $ZH_CLUSTER_ENV"
  unset _zh_dir _zh_cand _zh_found
elif [ "${ZH_ENV_OPTIONAL:-0}" = "1" ]; then
  echo "NOTE: no cluster.env found; continuing with built-in defaults (this step needs no SLURM"
  echo "      settings). Create one before submitting any GPU stage:"
  echo "          cp $_zh_dir/cluster.env.example $_zh_dir/cluster.env && \$EDITOR $_zh_dir/cluster.env"
  unset _zh_dir _zh_cand _zh_found
else
  echo "ERROR: cluster.env not found -- refusing to run with built-in defaults." >&2
  echo "  Looked at, in order:" >&2
  echo "    \$ZH_CLUSTER_ENV   ${ZH_CLUSTER_ENV:-<unset>}" >&2
  echo "    \$PWD              $PWD/cluster.env" >&2
  echo "    \$SLURM_SUBMIT_DIR ${SLURM_SUBMIT_DIR:+$SLURM_SUBMIT_DIR/cluster.env}${SLURM_SUBMIT_DIR:-<unset>}" >&2
  echo "    script directory  $_zh_dir/cluster.env" >&2
  echo "  cluster.env carries the conda env names, the system basename and the run lengths, so" >&2
  echo "  without it this would run against a different environment and still report success." >&2
  echo "  Create one from the template (gitignored: it is per-cluster and per-user):" >&2
  echo "      cp $_zh_dir/cluster.env.example $_zh_dir/cluster.env && \$EDITOR $_zh_dir/cluster.env" >&2
  echo "  Then submit with ./submit.sh, which locates it and exports ZH_CLUSTER_ENV to the job." >&2
  unset _zh_dir _zh_cand _zh_found
  return 1 2>/dev/null || exit 1
fi

#!/bin/bash
# Load the site settings (cluster.env) and, in a build directory, that build's sampling settings.
# *** SOURCE this file; do not execute it. ***
#
#     source "${SLURM_SUBMIT_DIR:-$PWD}/cluster_env.sh" || exit 1
#
# TWO layers, deliberately separate:
#
#   cluster.env   at the REPOSITORY ROOT. Site facts: account, partitions, resources, wall-times,
#                 conda env names. One per machine, shared by every build.
#   sampling.env  in each BUILD DIRECTORY, written by 01_build_system.sh. What changes the
#                 science: run lengths, replica count, the system basename. Per build, so
#                 re-running one system with different sampling cannot pick up another's.
#
# sampling.env is sourced AFTER cluster.env, so a build's choices win over site defaults.
#
# Every stage that reads a ZH_* setting REQUIRES cluster.env. Silently falling back to built-in
# defaults is worse than stopping: the job would activate a different conda env, or run under a
# different system basename, or for a different number of ns -- and still exit 0 looking like it
# worked. A wrong answer that reports success is the failure worth spending lines to prevent.
#
# The exception is 01_build_system.sh, which sets ZH_ENV_OPTIONAL=1: building needs no SLURM
# settings, so a system can be built before the cluster values are known.
#
# cluster.env is found by walking UP from the working directory, so it resolves from the bundle
# (../../cluster.env) and from a build directory (../../../cluster.env) without either needing to
# know how deep it sits. $ZH_CLUSTER_ENV, exported by submit.sh, always wins -- that is what makes
# a job use exactly what its submission used. Job scripts cannot locate it relative to themselves:
# sbatch copies the job script to a spool directory, so $0 there is not the repository.

_zh_up() {  # first directory at or above $1 containing cluster.env
  local d="$1" i
  for i in 1 2 3 4 5 6 7; do
    [ -f "$d/cluster.env" ] && { printf '%s\n' "$d/cluster.env"; return 0; }
    [ "$d" = "/" ] && break
    d="$(dirname "$d")"
  done
  return 1
}

_zh_found=""
for _zh_cand in "${ZH_CLUSTER_ENV:-}" "$(_zh_up "$PWD" || true)" \
                "${SLURM_SUBMIT_DIR:+$(_zh_up "$SLURM_SUBMIT_DIR" || true)}"; do
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
  echo "site settings:     $ZH_CLUSTER_ENV"
elif [ "${ZH_ENV_OPTIONAL:-0}" = "1" ]; then
  echo "NOTE: no cluster.env found; continuing with built-in defaults (this step needs no SLURM"
  echo "      settings). Create one at the repository root before submitting any batch stage:"
  echo "          cp cluster.env.example cluster.env && \$EDITOR cluster.env"
else
  echo "ERROR: cluster.env not found -- refusing to run with built-in defaults." >&2
  echo "  It belongs at the REPOSITORY ROOT and is found by walking up from the working" >&2
  echo "  directory. Looked from: $PWD${SLURM_SUBMIT_DIR:+ and $SLURM_SUBMIT_DIR}" >&2
  echo "  \$ZH_CLUSTER_ENV = ${ZH_CLUSTER_ENV:-<unset>}" >&2
  echo "  It carries the conda env names, partitions and wall-times, so without it this would" >&2
  echo "  run against a different environment and still report success. Create it (gitignored:" >&2
  echo "  per-cluster and per-user):" >&2
  echo "      cd <repository root> && cp cluster.env.example cluster.env && \$EDITOR cluster.env" >&2
  unset _zh_up _zh_cand _zh_found
  return 1 2>/dev/null || exit 1
fi

# Per-build sampling, when we are standing in a build directory.
if [ -f "$PWD/sampling.env" ]; then
  # shellcheck disable=SC1091
  . "$PWD/sampling.env"
  echo "sampling settings: $PWD/sampling.env"
fi
unset _zh_up _zh_cand _zh_found

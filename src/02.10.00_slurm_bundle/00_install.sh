#!/bin/bash
# Build the simulation conda env on the cluster. Run once on a login node.
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"

# Uses the repo's pinned spec. TODO(OQ-3): confirm the CUDA build matches the cluster driver.
conda env create -f ../../environment-cluster.yml || conda env update -f ../../environment-cluster.yml
conda activate zh853mor-sim

python - <<'PY'
import openmm, openmm.version
from openmm import Platform
print("OpenMM", openmm.version.version)
print("Platforms:", [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
PY
echo "Env ready. If CUDA is absent above, load the cluster CUDA module and reinstall openmm."

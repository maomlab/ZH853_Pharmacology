# ==============================================================================
# ZH853-MOR analysis. Environments (which interpreter each target needs):
#
#   LOCAL analysis env  -- EVERY target in this Makefile runs here.
#                          Create with `make env-local`, then `conda activate zh853mor-local`.
#                          (`molstar-render` additionally needs Node.js >= 18.)
#
#   CLUSTER envs (zh853mor-prep / -sim / -plumed) -- used ONLY by the SLURM bundle
#                          in src/02.10.00_slurm_bundle/, run ON the cluster, NOT via
#                          make. See that directory's README. `make env-cluster`
#                          just creates those envs if you are on the cluster.
#
# Env specs are named environment_<env name>.yml, one per conda env.
#
# Dependencies are encoded as prerequisites, so `make <target>` first runs whatever
# it needs (e.g. `make interactions` runs `fetch`; `make molstar-render` runs
# `prep-complex-split`). Prereqs are cheap/idempotent (fetch skips existing files).
# Run `make help` for the grouped target list.
# ==============================================================================

.PHONY: help \
        env-local env-cluster \
        lint format typecheck test check \
        fetch \
        qc interactions mutations analogs design interaction-map depictions analysis \
        prep-complex-split prep-receptor-protonate prep-receptor-rebuild \
        prep-ZH853-protonate prep-receptor-orient prep-analogs-pose \
        prep-ligand-parameterize prep \
        membrane-plot \
        molstar-render figures manuscript clean-intermediate

help:  ## Show this grouped target list
	@awk 'BEGIN {FS = ":.*?## "} \
		/^## / {printf "\n\033[1m%s\033[0m\n", substr($$0, 4)} \
		/^[a-zA-Z0-9_-]+:.*?## / {printf "  \033[36m%-18s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

## Environments
env-local:  ## Create the LOCAL analysis env zh853mor-local (this Makefile runs here)
	conda env create -f environment_zh853mor-local.yml || conda env update -f environment_zh853mor-local.yml

env-cluster:  ## Create ALL cluster envs: zh853mor-prep, -sim, -plumed (on the cluster)
	conda env create -f environment_zh853mor-prep.yml   || conda env update -f environment_zh853mor-prep.yml
	conda env create -f environment_zh853mor-sim.yml    || conda env update -f environment_zh853mor-sim.yml
	conda env create -f environment_zh853mor-plumed.yml || conda env update -f environment_zh853mor-plumed.yml

## Development & CI  [local env]
lint:  ## Ruff lint (package + tests)
	ruff check src/zh853mor tests

format:  ## Ruff auto-format
	ruff format src tests
	ruff check --fix src/zh853mor tests

typecheck:  ## Mypy type-check
	mypy

test:  ## Run unit tests
	pytest

check: lint typecheck test  ## Lint + typecheck + test (the CI gate)

## Data acquisition  [local env]
fetch:  ## Download comparator PDBs -> data/comparators/ (prereq of: interactions, mutations)
	python src/01.01.00_fetch_comparators.py

## Static analysis - Objectives 1-3  [local env]
qc:  ## Structure QC report -> product/
	python src/01.02.00_qc_structure.py

interactions: fetch  ## Comparative interaction fingerprints + heatmap -> product/ (Obj 1)
	python src/03.01.00_interaction_fingerprints.py

mutations: fetch  ## Ranked ZH853-selective mutation panel -> product/ (Obj 2)
	python src/03.02.00_mutation_panel.py

analogs:  ## Analog physicochemical property panel -> product/ (Obj 3)
	python src/05.01.00_analog_properties.py

design:  ## Structure-guided modification design -> product/ (Obj 3)
	python src/05.02.00_design_modifications.py

interaction-map:  ## PoseView-style 2D interaction map -> manuscript fig5 (Obj 1)
	python src/03.03.00_interaction_map.py

depictions:  ## 2D vector ligand depictions -> manuscript fig4 (Obj 3)
	python src/05.03.00_ligand_depictions.py

analysis: qc interactions mutations analogs design  ## Run the full static-analysis pipeline

## MD system prep - Phase 2  [local env, EXCEPT prep-ligand-parameterize; feeds the SLURM bundle]
# Targets are prep-<object>-<action> and are listed in script order, so `make help` reads as the
# running order. The receptor is oriented (02.05.00) AFTER the ligand is prepared (02.04.00)
# because that is the script numbering; the two are independent.
prep-complex-split:  ## 02.01.00  Assess prep needs + split the complex into components -> product/, intermediate/
	python src/02.01.00_assess_and_split.py

prep-receptor-protonate:  ## 02.02.00  PROPKA protonation states for the receptor -> product/
	python src/02.02.00_protonation.py

prep-receptor-rebuild:  ## 02.03.00  Rebuild receptor sidechains + caps (PDBFixer) -> intermediate/
	python src/02.03.00_prepare_receptor.py

prep-ZH853-protonate:  ## 02.04.00  Protonated ZH853 (+1) + parameterization inputs -> intermediate/
	python src/02.04.00_ligand_prep.py

prep-receptor-orient: prep-receptor-rebuild  ## 02.05.00  Superpose onto OPM so the membrane normal is z (PACKMOL-Memgen --preoriented)
	python src/02.05.00_orient_receptor.py

# No prerequisites ON PURPOSE. Its only inputs are intermediate/02.05.00_oriented/complex_oriented.pdb
# and rdkit, and the script checks for that itself. Depending on prep-receptor-orient would drag in
# prep-receptor-rebuild, which needs pdbfixer/openmm -- absent from the cluster prep env by design --
# so this target would be unrunnable on the cluster for a step that has no such requirement.
# `make prep` still runs the whole chain in order.
prep-analogs-pose:  ## 02.07.00  ZH850/ZH831/ZH809 poses by scaffold transfer from ZH853 (also runs in zh853mor-prep on the cluster)
	python src/02.07.00_analog_poses.py

prep-ligand-parameterize:  ## 02.08.00  GAFF2/AM1-BCC parameters for every ligand -> intermediate/ [needs AmberTools: zh853mor-prep, NOT zh853mor-local]
	bash src/02.08.00_ligand_parameterize.sh

# prep-ligand-parameterize is deliberately NOT here: it needs AmberTools (zh853mor-prep), and
# every other target in this Makefile runs in zh853mor-local. Run it separately, in that env.
prep: prep-complex-split prep-receptor-protonate prep-receptor-rebuild prep-ZH853-protonate \
      prep-receptor-orient prep-analogs-pose  ## Run the full Phase-2 local prep, in script order

## Figures & manuscript  [LOCAL env (zh853mor-local); molstar-render also needs Node.js >= 18]
membrane-plot: prep-receptor-orient  ## 03.04.00  Membrane-placement determination plot -> product/ (manuscript panel B)
	python src/03.04.00_membrane_placement.py

molstar-render: prep-complex-split prep-receptor-orient  ## Headless MolStar 3D renders -> manuscript (complex, pocket, membrane; needs Node.js)
	cd src/03.10.00_molstar_render && npm install --silent \
	  && python build_overview_camera.py && python build_pocket_mvs.py && python build_membrane_mvs.py \
	  && node render.js && python trim_figures.py

figures: interactions analogs design depictions interaction-map membrane-plot molstar-render  ## Regenerate + stage ALL manuscript figures
	cp $$(ls -t product/03.01.00_fingerprint_heatmap_*.png  | head -1) product/manuscript/figures/fig1_interaction_heatmap.png
	cp $$(ls -t product/05.01.00_analog_property_space_*.png | head -1) product/manuscript/figures/fig2_property_space.png
	cp $$(ls -t product/05.02.00_design_property_shifts_*.png | head -1) product/manuscript/figures/fig3_design_shifts.png

manuscript:  ## Compile the LaTeX -> manuscript.pdf (run `make figures` first if figures changed)
	cd product/manuscript && tectonic manuscript.tex

## Housekeeping
clean-intermediate:  ## Remove cached intermediate results
	rm -rf intermediate/*

# ==============================================================================
# ZH853-MOR analysis. Environments (which interpreter each target needs):
#
#   LOCAL analysis env  -- EVERY target in this Makefile runs here.
#                          Create with `make env-local`, then `conda activate zh853mor-local`.
#
#   NODE (>= 18)        -- the two headless-MolStar render stages, `molstar-render` (03.10.00) and
#                          `movie-molstar` (02.12.00), need it plus an installed node_modules.
#                          Locally both targets run `npm install` themselves. On the CLUSTER node
#                          ships inside zh853mor-prep and the install must be done on a LOGIN
#                          node: `make env-cluster` does it, or `make env-node` on its own.
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
        env-local env-cluster env-node \
        lint format typecheck test check \
        fetch \
        qc interactions mutations analogs design interaction-map depictions analysis \
        prep-complex-split prep-receptor-protonate prep-receptor-rebuild \
        prep-ZH853-protonate prep-receptor-orient prep-analogs-pose \
        prep-ligand-parameterize prep \
        sim-reduce sim-aggregate sim-figures simulation-analysis \
        movie-export movie-dashboard movie-molstar movies \
        landscape-fit landscape-figures landscape \
        membrane-plot \
        molstar-render figures manuscript clean-intermediate

help:  ## Show this grouped target list
	@awk 'BEGIN {FS = ":.*?## "} \
		/^## / {printf "\n\033[1m%s\033[0m\n", substr($$0, 4)} \
		/^[a-zA-Z0-9_-]+:.*?## / {printf "  \033[36m%-18s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

## Environments
env-local:  ## Create the LOCAL analysis env zh853mor-local (this Makefile runs here)
	conda env create -f environment_zh853mor-local.yml || conda env update -f environment_zh853mor-local.yml

env-cluster:  ## Create ALL cluster envs + the Node deps (RUN ON A LOGIN NODE)
	conda env create -f environment_zh853mor-prep.yml   || conda env update -f environment_zh853mor-prep.yml
	conda env create -f environment_zh853mor-sim.yml    || conda env update -f environment_zh853mor-sim.yml
	conda env create -f environment_zh853mor-plumed.yml || conda env update -f environment_zh853mor-plumed.yml
	@echo ""
	@echo "--- Node deps for the headless MolStar movie pass (02.12.00) ---"
	@cd src/02.12.00_render_trajectory_movies \
	  && conda run -n zh853mor-prep npm install --silent \
	  && echo "  installed: the MolStar movie pass will run on this cluster" \
	  || { echo "  SKIPPED: npm install failed."; \
	       echo "  Usually no outbound network -- run 'make env-cluster' on a LOGIN node."; \
	       echo "  Without it submit_render.sh skips the MolStar pass with a message; the"; \
	       echo "  dashboard movies and every other stage are unaffected."; }

# The two headless-MolStar renderers each keep their own node_modules, because puppeteer downloads
# its own Chromium (~150 MB apiece). The LOCAL targets install on demand -- `movie-molstar` and
# `molstar-render` both run `npm install` themselves -- so this target exists for the CLUSTER,
# where the install has to happen on a LOGIN node: a compute node usually has no outbound network,
# and a job that discovers that has already spent its queue time.
NODE_STAGES := src/02.12.00_render_trajectory_movies src/03.10.00_molstar_render

env-node:  ## Install the headless-MolStar Node deps for both render stages (LOGIN NODE)
	@command -v node >/dev/null 2>&1 || { \
	  echo "ERROR: node is not on PATH."; \
	  echo "  Cluster: it ships with zh853mor-prep -- 'conda activate zh853mor-prep', then retry."; \
	  echo "  Local:   install Node.js >= 18."; exit 1; }
	@node --version | sed 's/^/node /'
	@for d in $(NODE_STAGES); do \
	  echo "npm install in $$d"; \
	  (cd $$d && npm install --silent) \
	    || echo "  WARNING: failed in $$d (on a cluster, usually no outbound network)"; \
	done

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

## Production MD analysis - Phase 4  [local env; step 01 runs ON THE CLUSTER]
# 01_reduce_trajectory.py is deliberately NOT a make target of its own beyond `sim-reduce`: it
# reads the trajectories, which stay on the cluster (~100 GB for the panel), and runs there as a
# job array (src/02.11.00_analyze_simulations/submit_reduce.sh). Everything below it works on the
# few hundred kB per replica that reduction writes, so it runs locally.
sim-reduce:  ## 02.11.00  Reduce production replicas -> intermediate/ [CLUSTER: needs the DCDs + zh853mor-prep]
	python src/02.11.00_analyze_simulations/01_reduce_trajectory.py --all

sim-aggregate:  ## 02.11.00  QC verdicts + occupancy tables + report -> product/
	python src/02.11.00_analyze_simulations/02_aggregate.py

sim-figures:  ## 02.11.00  QC dashboard, convergence, occupancy, pocket dynamics -> product/
	python src/02.11.00_analyze_simulations/03_figures.py

simulation-analysis: sim-aggregate sim-figures  ## Full local MD analysis (after sim-reduce on the cluster)

## Trajectory movies - Phase 4  [local env; step 01 runs ON THE CLUSTER; movie-molstar needs Node.js >= 18]
# Same cluster/local split as 02.11.00 and for the same reason: movie-export reads the DCDs, which
# stay on the cluster, and writes the few tens of MB per replica that the two renders work from.
movie-export:  ## 02.12.00  Export viewable movie trajectories -> intermediate/ [CLUSTER: needs the DCDs + zh853mor-prep]
	python src/02.12.00_render_trajectory_movies/01_export_movie_trajectory.py --all

movie-dashboard:  ## 02.12.00  Diagnostic movies: structure + its own time series on one clock -> product/
	python src/02.12.00_render_trajectory_movies/02_render_dashboard.py

movie-molstar:  ## 02.12.00  Cartoon MolStar movies -> product/ (needs Node.js)
	cd src/02.12.00_render_trajectory_movies && npm install --silent && node 03_render_molstar.js

# Rendering is serial WITHIN a replica (matplotlib -> ffmpeg frame by frame; MolStar one frame at
# a time in software WebGL), so a panel takes N times one replica. These targets render the whole
# panel serially, which is right for one or two replicas; for the whole panel use the job array,
# src/02.12.00_render_trajectory_movies/submit_render.sh, which runs one replica per task.
movies: movie-dashboard movie-molstar  ## Both movie renders, SERIALLY (cluster: use submit_render.sh)

## Conformational landscapes - Phase 4  [local env; needs 02.11.00's reduced replicas]
# Entirely local: the per-frame feature matrices these fit on were written by 02.11.00's
# reduction, which is the only stage that reads the trajectories. SPACE selects the coordinates
# (receptor | ligand | pair); pass through for the others, e.g. `make landscape SPACE=ligand`.
SPACE ?= receptor

landscape-fit:  ## 02.13.00  Fit the shared tICA basis + project every replica -> intermediate/
	python src/02.13.00_map_conformational_landscapes/01_fit_landscape.py --space $(SPACE) --scan

landscape-figures:  ## 02.13.00  Landscapes, difference maps and diagnostics -> product/
	python src/02.13.00_map_conformational_landscapes/02_figures.py --space $(SPACE)

landscape: landscape-fit landscape-figures  ## Both landscape steps (SPACE=receptor|ligand|pair)

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

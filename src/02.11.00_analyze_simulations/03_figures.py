#!/usr/bin/env python3
"""Figures for the production MD: QC dashboard, convergence, contacts, pocket dynamics.

**Runs locally** (`zh853mor-local`) on the reduced replicas written by `01_reduce_trajectory.py`.
Four figures, one per question:

``qc_dashboard``       is the simulation physically sound and structurally stable?
``convergence``        is it sampled well enough to average? (blocking, running means, cosine)
``contact_occupancy``  which pocket contacts persist, per system, with replicate spread?
``pocket_dynamics``    the Objective-1/2 and D-11 comparisons: anchors, activation, Na+ site

Every panel that needs data not present is skipped with a note rather than drawn empty, so the
figures are usable while only part of the panel of systems has run.

Run: ``python src/02.11.00_analyze_simulations/03_figures.py``  (or ``make sim-figures``).
"""

from __future__ import annotations

import argparse
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from zh853mor import convergence as cv  # noqa: E402
from zh853mor import md, paths, structure  # noqa: E402

ANALYSIS_ROOT = paths.INTERMEDIATE / "02.11.00_analyze_simulations"
PREFIX = "02.11.00"
NA_SITE_CUT = 3.2  # A; matches 02_aggregate.py


def colors(systems: list[str]) -> dict[str, tuple]:
    """One colour per system, fixed by position so every figure agrees."""
    cmap = plt.get_cmap("tab10")
    return {s: cmap(i % 10) for i, s in enumerate(systems)}


def empty(ax, message: str) -> None:
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=7, style="italic",
            transform=ax.transAxes, wrap=True)
    ax.set_xticks([])
    ax.set_yticks([])


def plot_series(ax, systems, col, name, ylabel, title):
    """Every replica as a thin line, in the system's colour; nothing if the series is absent."""
    drawn = 0
    labelled: set[str] = set()   # one legend entry per SYSTEM, not per replica
    for s, reps in systems.items():
        for r in reps:
            try:
                y = r.series(name, equilibrated=False)
                t = r.series("time_ns", equilibrated=False)
            except KeyError:
                continue
            if not np.isfinite(y).any():
                continue
            ax.plot(t, y, lw=0.6, alpha=0.8, color=col[s],
                    label=None if s in labelled else s)
            labelled.add(s)
            # Mark where the relaxation was judged to end; means are taken after this.
            if 0 < r.t0 < len(t):
                ax.axvline(t[r.t0], color=col[s], lw=0.5, ls=":", alpha=0.5)
            drawn += 1
    if not drawn:
        empty(ax, f"no {name}")
        return
    ax.set_xlabel("time (ns)", fontsize=7)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.set_title(title, fontsize=8)
    ax.tick_params(labelsize=6)


def qc_dashboard(systems, col, out: Path) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(13, 6.5))
    plot_series(axes[0, 0], systems, col, "rmsd_ca_tm", "RMSD (A)",
                "TM Ca RMSD vs staged receptor")
    plot_series(axes[0, 1], systems, col, "lig_rmsd_pose", "RMSD (A)",
                "ligand vs deposited pose (receptor-aligned)")

    ax = axes[0, 2]
    drawn = False
    stale = []
    for s, reps in systems.items():
        stacks, resids = [], None
        for r in reps:
            arr = r.arrays()
            if "rmsf" not in arr or not arr["rmsf"].size:
                continue
            # RMSF has one value per Ca, which is NOT one per residue: the ACE/NME caps have
            # none. `ca_resids` is that axis. Files written before it existed pair a 281-value
            # RMSF with 283 residue ids, so say what to do rather than raising a shape error.
            x = arr.get("ca_resids")
            if x is None or x.size != arr["rmsf"].size:
                stale.append(f"{s}/{r.replica}")
                continue
            stacks.append(arr["rmsf"])
            resids = x
        if stacks and resids is not None:
            m = np.mean(stacks, axis=0)
            sd = np.std(stacks, axis=0)
            ax.plot(resids, m, lw=0.8, color=col[s], label=s)
            ax.fill_between(resids, m - sd, m + sd, color=col[s], alpha=0.2, lw=0)
            drawn = True
    if drawn:
        for resid in md.ANCHORS:
            ax.axvline(resid, color="0.6", lw=0.4, ls=":")
        ax.set_xlabel("residue (human OPRM1)", fontsize=7)
        ax.set_ylabel("RMSF (A)", fontsize=7)
        ax.set_title("Ca RMSF (post-equilibration; dotted = anchors)", fontsize=8)
        ax.tick_params(labelsize=6)
    else:
        empty(ax, "no RMSF" if not stale else
              "RMSF needs re-reduction:\nthese files predate `ca_resids`\n"
              "(rerun 01_reduce_trajectory.py --force)")
    if stale:
        print(f"  note: RMSF skipped for {len(stale)} replica(s) written before the Ca residue "
              f"axis was stored ({', '.join(stale[:3])}"
              f"{', ...' if len(stale) > 3 else ''}); rerun the reduction with --force.")

    plot_series(axes[1, 0], systems, col, "apl_net", "A^2/lipid",
                "area per lipid (protein-corrected)")
    plot_series(axes[1, 1], systems, col, "thickness", "A", "bilayer thickness (P-P)")
    plot_series(axes[1, 2], systems, col, "log_density", "g/mL", "density")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=min(len(labels), 6), fontsize=7,
                   frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("Production MD quality control (dotted vertical = end of equilibration)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def convergence_figure(systems, col, out: Path) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))

    ax = axes[0]
    drawn = False
    for s, reps in systems.items():
        for r in reps:
            try:
                y = r.series("rmsd_ca_tm")
            except KeyError:
                continue
            y = y[np.isfinite(y)]
            if y.size < 20:
                continue
            sizes, sems = cv.block_sem_curve(y)
            ax.plot(sizes, sems, lw=0.8, color=col[s], alpha=0.8)
            drawn = True
    if drawn:
        ax.set_xscale("log")
        ax.set_xlabel("block size (frames)", fontsize=7)
        ax.set_ylabel("SEM of the block means (A)", fontsize=7)
        ax.set_title("Blocking: TM Ca RMSD\n(plateau = the honest error bar)", fontsize=8)
        ax.tick_params(labelsize=6)
    else:
        empty(ax, "no TM RMSD series")

    ax = axes[1]
    drawn = False
    for s, reps in systems.items():
        for r in reps:
            try:
                y = r.series("rmsd_ca_tm")
                t = r.series("time_ns")
            except KeyError:
                continue
            mask = np.isfinite(y)
            y, t = y[mask], t[mask]
            if y.size < 20:
                continue
            running = np.cumsum(y) / np.arange(1, y.size + 1)
            ax.plot(t, running, lw=0.8, color=col[s], alpha=0.8)
            drawn = True
    if drawn:
        ax.set_xlabel("time (ns)", fontsize=7)
        ax.set_ylabel("running mean (A)", fontsize=7)
        ax.set_title("Running mean of the TM Ca RMSD\n(flat and superposed = converged)",
                     fontsize=8)
        ax.tick_params(labelsize=6)
    else:
        empty(ax, "no TM RMSD series")

    ax = axes[2]
    names, values, errors = [], [], []
    for s, reps in systems.items():
        cos = [r.summary.get("pca", {}).get("cosine_content", [np.nan])[0] for r in reps]
        cos = [float(c) for c in cos if np.isfinite(float(c))]
        if cos:
            names.append(s)
            values.append(float(np.mean(cos)))
            errors.append(float(np.std(cos, ddof=1) / np.sqrt(len(cos))) if len(cos) > 1 else 0.0)
    if names:
        ax.bar(range(len(names)), values, yerr=errors, capsize=3,
               color=[col[n] for n in names])
        ax.axhline(0.5, color="crimson", lw=1, ls="--")
        ax.text(0.02, 0.52, "above: PC1 is diffusion (Hess 2002)", fontsize=6, color="crimson",
                transform=ax.get_yaxis_transform())
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=30, ha="right", fontsize=6)
        ax.set_ylim(0, 1)
        ax.set_ylabel("PC1 cosine content", fontsize=7)
        ax.set_title("Is the leading motion real?", fontsize=8)
        ax.tick_params(labelsize=6)
    else:
        empty(ax, "no PCA")

    fig.suptitle("Convergence diagnostics (each line is one replica)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def occupancy_figure(systems, col, out: Path) -> Path | None:
    """Heatmap of pocket-residue occupancy + anchor bars with replicate error bars."""
    holo = {s: reps for s, reps in systems.items()
            if any(r.summary.get("n_ligand_atoms") for r in reps)}
    if not holo:
        return None
    resids: np.ndarray | None = None
    mat, labels = [], []
    for s, reps in holo.items():
        stacks = []
        for r in reps:
            arr = r.arrays()
            if "occupancy" in arr and arr["occupancy"].size:
                stacks.append(arr["occupancy"])
                resids = arr["resids"]
        if stacks:
            mat.append(np.mean(stacks, axis=0))
            labels.append(s)
    if not mat or resids is None:
        return None
    occ = np.array(mat)

    # Show the pocket shell (the static Objective-1 residue set) plus anything the MD finds
    # occupied above 20% -- a contact the 3.5 A structure does not have is a result, not noise.
    keep = sorted({int(r) for r in resids if int(r) in md.KEY_RESIDUES}
                  | {int(r) for i, r in enumerate(resids) if occ[:, i].max() >= 0.2})
    idx = [int(np.flatnonzero(resids == r)[0]) for r in keep]

    fig, axes = plt.subplots(1, 2, figsize=(max(8.0, 0.22 * len(keep) + 4), 5),
                             gridspec_kw={"width_ratios": [3, 1.4]})
    ax = axes[0]
    im = ax.imshow(occ[:, idx], aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xticks(range(len(keep)))
    ax.set_xticklabels([f"{structure.bw(r) if structure.bw(r) != '-' else ''} {r}".strip()
                        for r in keep], rotation=90, fontsize=5)
    ax.set_title("Contact occupancy (fraction of frames within "
                 f"{md.CONTACT_CUT} A), mean over replicas", fontsize=8)
    cb = fig.colorbar(im, ax=ax, shrink=0.6)
    cb.set_label("occupancy", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    ax = axes[1]
    anchors = list(md.ANCHORS)
    width = 0.8 / max(len(labels), 1)
    for k, s in enumerate(labels):
        vals, errs = [], []
        for a in anchors:
            per_rep = [float(r.summary.get("occupancy", {}).get(str(a), np.nan))
                       for r in holo[s]]
            per_rep = [v for v in per_rep if np.isfinite(v)]
            vals.append(float(np.mean(per_rep)) if per_rep else np.nan)
            errs.append(float(np.std(per_rep, ddof=1) / np.sqrt(len(per_rep)))
                        if len(per_rep) > 1 else 0.0)
        ax.barh(np.arange(len(anchors)) + k * width, vals, height=width, xerr=errs,
                color=col[s], label=s, capsize=2, error_kw={"lw": 0.7})
    ax.set_yticks(np.arange(len(anchors)) + 0.4 - width / 2)
    ax.set_yticklabels([f"{md.ANCHORS[a]} ({a})" for a in anchors], fontsize=6)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("occupancy", fontsize=7)
    ax.set_title("Anchors, +- SEM over replicas", fontsize=8)
    ax.legend(fontsize=6, frameon=False)
    ax.tick_params(labelsize=6)

    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def pocket_dynamics(systems, col, out: Path) -> Path:
    """The comparisons the objectives ask for: pose, activation, and the D2.50 sodium site."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))

    ax = axes[0]
    drawn = False
    for s, reps in systems.items():
        pooled = []
        for r in reps:
            try:
                y = r.series("lig_rmsd_pose")
            except KeyError:
                continue
            pooled.append(y[np.isfinite(y)])
        if pooled and sum(p.size for p in pooled):
            ax.hist(np.concatenate(pooled), bins=40, histtype="step", density=True,
                    color=col[s], label=s, lw=1.0)
            drawn = True
    if drawn:
        ax.set_xlabel("ligand RMSD vs deposited pose (A)", fontsize=7)
        ax.set_ylabel("density", fontsize=7)
        ax.set_title("Pose retention", fontsize=8)
        ax.legend(fontsize=6, frameon=False)
        ax.tick_params(labelsize=6)
    else:
        empty(ax, "no ligand system reduced yet")

    ax = axes[1]
    names, means, errs = [], [], []
    for s, reps in systems.items():
        vals = [r.mean_of("tm3_tm6_ca") for r in reps]
        vals = [v for v in vals if np.isfinite(v)]
        if vals:
            names.append(s)
            means.append(float(np.mean(vals)))
            errs.append(float(np.std(vals, ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else 0.0)
    if names:
        ax.bar(range(len(names)), means, yerr=errs, capsize=3, color=[col[n] for n in names])
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=30, ha="right", fontsize=6)
        ax.set_ylabel("R3.50-T6.34 Ca distance (A)", fontsize=7)
        ax.set_title("TM6 opening: does the ligand hold\nthe active state without a transducer?",
                     fontsize=8)
        ax.tick_params(labelsize=6)
    else:
        empty(ax, "no activation ruler")

    ax = axes[2]
    drawn = False
    for s, reps in systems.items():
        pooled = []
        for r in reps:
            try:
                d = r.series("na_d250_dist")
            except KeyError:
                continue
            pooled.append(d[np.isfinite(d)])
        if pooled and sum(p.size for p in pooled):
            ax.hist(np.concatenate(pooled), bins=40, histtype="step", density=True,
                    color=col[s], label=s, lw=1.0,
                    ls="--" if s.endswith("ASH") else "-")
            drawn = True
    if drawn:
        ax.axvline(NA_SITE_CUT, color="crimson", lw=1, ls=":")
        ax.set_xlabel("Na+ to nearest D2.50 carboxylate O (A)", fontsize=7)
        ax.set_ylabel("density", fontsize=7)
        ax.set_title("D2.50 sodium site (D-11)\ndashed = ASH (protonated)", fontsize=8)
        ax.legend(fontsize=6, frameon=False)
        ax.tick_params(labelsize=6)
    else:
        empty(ax, "no sodium series")

    fig.suptitle("Pocket and state observables (pooled over replicas)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--analysis-root", type=Path, default=ANALYSIS_ROOT)
    ap.add_argument("--date", default=f"{date.today():%Y%m%d}")
    args = ap.parse_args()

    results = md.load_replicas(args.analysis_root)
    if not results:
        raise SystemExit(f"ERROR: no reduced replicas under {args.analysis_root}. "
                         "Run 01_reduce_trajectory.py on the cluster first.")
    systems = md.group_by_system(results)
    col = colors(list(systems))
    paths.ensure_dir(paths.PRODUCT)

    written = [
        qc_dashboard(systems, col, paths.PRODUCT / f"{PREFIX}_qc_dashboard_{args.date}.png"),
        convergence_figure(systems, col, paths.PRODUCT / f"{PREFIX}_convergence_{args.date}.png"),
        occupancy_figure(systems, col,
                         paths.PRODUCT / f"{PREFIX}_contact_occupancy_{args.date}.png"),
        pocket_dynamics(systems, col, paths.PRODUCT / f"{PREFIX}_pocket_dynamics_{args.date}.png"),
    ]
    for path in written:
        print(f"  {path}" if path else "  (contact occupancy skipped: no ligand system)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

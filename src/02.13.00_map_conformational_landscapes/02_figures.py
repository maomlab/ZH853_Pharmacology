#!/usr/bin/env python3
"""Free-energy landscapes on the shared 2D basis, per system, with difference maps and diagnostics.

**Runs locally** (`zh853mor-local`) on what `01_fit_landscape.py` wrote. Three figures:

``landscapes``   one -kT ln P heatmap per system, on IDENTICAL bins and one colour scale
``differences``  each system minus a reference, so the comparison is a picture rather than
                 two pictures the eye has to subtract
``diagnostics``  the four things that decide whether the landscapes mean anything: implied
                 timescale against lag, replica-to-replica agreement, cosine content, and what
                 the axes are made of

Every panel is binned on the same grid and clipped to the same energy range, because two
landscapes drawn on different bins or different colour limits are not comparable however alike
they look -- which is the whole purpose of the exercise.

    python 02_figures.py --space receptor
    python 02_figures.py --space ligand --reference ZH853_ASP --fmax 5
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm  # noqa: E402

from zh853mor import landscape as ls  # noqa: E402
from zh853mor import paths  # noqa: E402

LAND_ROOT = paths.INTERMEDIATE / "02.13.00_map_conformational_landscapes"
PREFIX = "02.13.00"

FMAX = 4.0             # kcal/mol; ~6.5 kT at 310 K, past which the density is a handful of frames
MIN_COUNT = 5          # frames per bin below which -kT ln P is one sample pretending to be a wall
COSINE_FLAG = 0.5      # D-20's threshold, applied to a tIC projection

SURFACE = "#fcfcfb"
UNSAMPLED = "#eae8e4"   # warm grey: never produced by the blue ramp, so "no data" cannot be read
                        # as "high free energy"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#dcdcd8"

# Sequential: ONE hue, the validated blue ramp (100 -> 700), used DARK = LOW free energy.
# That inverts the usual "lightest means near zero", deliberately: here near-zero free energy is
# the most-populated and best-determined region, and the pale end is the sparse rim where the
# estimate rests on a few frames. Putting the ink on the rim would weight the noise.
FES_CMAP = LinearSegmentedColormap.from_list(
    "fes", ["#0d366b", "#184f95", "#256abf", "#3987e5", "#6da7ec", "#9ec5f4", "#cde2fb"])
# Diverging: blue <-> red through a NEUTRAL grey, so "no difference" reads as nothing.
DIFF_CMAP = LinearSegmentedColormap.from_list(
    "delta", ["#0d366b", "#2a78d6", "#a8c8ee", "#f0efec", "#eda79a", "#d6452a", "#6b170d"])
# Replica identity, categorical slots 1-4 (validated: CVD dE 9.1, normal dE 22.9). The contrast
# WARN on slots 3-4 is discharged by the legend that always accompanies them.
REPLICA_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]


def load(space: str, root: Path) -> tuple[dict, dict[str, list[tuple[str, np.ndarray]]]]:
    d = root / space
    model_path, proj_path = d / "model.json", d / "projections.npz"
    if not (model_path.exists() and proj_path.exists()):
        raise SystemExit(f"ERROR: {d} has no fitted landscape. Run:\n"
                         f"    python 01_fit_landscape.py --space {space}")
    model = json.loads(model_path.read_text())
    proj: dict[str, list[tuple[str, np.ndarray]]] = {}
    with np.load(proj_path) as data:
        for key in data.files:
            if not key.startswith("proj__"):
                continue
            _, system, replica = key.split("__", 2)
            proj.setdefault(system, []).append((replica, np.asarray(data[key], dtype=float)))
    for reps in proj.values():
        reps.sort()
    return model, proj


def style(ax, xlabel: str = "", ylabel: str = "", title: str = "") -> None:
    ax.set_facecolor(UNSAMPLED)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=7, color=MUTED)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=7, color=MUTED)
    if title:
        ax.set_title(title, fontsize=8, loc="left", color=INK)
    ax.tick_params(labelsize=6, length=2, colors=MUTED)
    for s in ax.spines.values():
        s.set_color(GRID)


def grid_shape(n: int) -> tuple[int, int]:
    cols = min(4, max(1, int(np.ceil(np.sqrt(n)))))
    return int(np.ceil(n / cols)), cols


def edges(proj: dict[str, list[tuple[str, np.ndarray]]], bins: int):
    """One grid for every system -- the precondition for the panels being comparable at all."""
    allp = [p for reps in proj.values() for _, p in reps]
    xe = np.linspace(*ls.common_extent([p[:, 0] for p in allp]), bins + 1)
    ye = np.linspace(*ls.common_extent([p[:, 1] for p in allp]), bins + 1)
    return xe, ye


def surfaces(proj, xe, ye, min_count: int) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    out = {}
    for system, reps in proj.items():
        p = np.vstack([q for _, q in reps])
        out[system] = ls.free_energy(p[:, 0], p[:, 1], xe, ye, min_count=min_count)
    return out


def axis_labels(model: dict) -> tuple[str, str]:
    """Axis names, with units where the feature's prefix makes the unit certain."""
    names = model["components"][:2]
    if model["space"] == "pair":
        # `ca:` and `lig:` are distances by construction; a bare series name could be anything
        # the reduction stored (an RMSD, a thickness, a radius of gyration), so it is left alone
        # rather than labelled with a unit that might be wrong.
        return tuple(f"{n} (Å)" if n.startswith(("ca:", "lig:", "activation:")) else n
                     for n in names)  # type: ignore[return-value]
    return f"{names[0]}  (kinetic map)", f"{names[1]}  (kinetic map)"


def header(fig, model: dict, what: str, note: str = "", colorbar: bool = True) -> None:
    """Title block placed in INCHES from the top, not in figure fractions.

    A fraction that clears the title on a two-row grid runs straight through it on a four-row
    one, because the same 0.045 is a different number of points at a different figure height.
    """
    h = fig.get_figheight()
    lag = model["lag"]
    fig.text(0.012, 1 - 0.26 / h,
             f"{what} — {model['space']} space, basis fitted on {model['fit_on']}",
             fontsize=11, ha="left", va="top", color=INK)
    bits = [f"lag {lag['ns']} ns", f"{len(model['systems'])} systems",
            f"-kT ln P at {ls.TEMPERATURE_K:.0f} K, bins shared across panels"]
    tica = model.get("tica")
    if tica:
        bits.insert(1, f"{tica['n_features']} features → rank {tica['rank']}")
    if note:
        bits.append(note)
    fig.text(0.012, 1 - 0.52 / h, " | ".join(bits), fontsize=7, ha="left", va="top", color=MUTED)
    fig.subplots_adjust(top=1 - 0.80 / h, bottom=0.78 / h, left=0.07,
                        right=0.88 if colorbar else 0.98)


def footnote(fig, text: str) -> None:
    """Caveat text, wrapped to the figure's actual width.

    matplotlib will not wrap a `fig.text`, and the panel grid is sized from the system count, so
    a one-line caveat that fits a four-panel figure runs off the edge of a two-panel one -- and
    the half that gets cut is the half that says what the picture does not mean.
    """
    import textwrap
    width = max(60, int(fig.get_figwidth() * 72 / 3.9))     # ~3.9 pt per char at fontsize 6.5
    lines = textwrap.wrap(text, width=width)
    fig.text(0.012, 0.16 / fig.get_figheight(), "\n".join(lines), fontsize=6.5, color=MUTED,
             ha="left", va="bottom", linespacing=1.5)


def auto_bins(proj, requested: int | None) -> int:
    """Bins per axis, scaled to the sample size unless the caller fixed it.

    A grid fine enough for 15 000 frames leaves a 900-frame replica with a handful of frames in
    every occupied bin and a landscape made of single-sample specks. Aiming at roughly 8 frames
    per occupied bin, and assuming the occupied region covers about a quarter of the bounding
    box, gives bins ~ sqrt(n / 2); the clamp keeps it legible at both extremes.
    """
    if requested:
        return requested
    per_system = [sum(len(p) for _, p in reps) for reps in proj.values()]
    n = float(np.median(per_system)) if per_system else 0.0
    return int(np.clip(round(np.sqrt(max(n, 1.0) / 2.0)), 16, 60))


def hdr_contour(ax, p: np.ndarray, xe, ye, fraction: float, color: str, ls_: str):
    """Contour enclosing `fraction` of one replica's frames -- a highest-density region.

    Used instead of a free-energy contour for the per-replica panel: a single replica has a few
    hundred post-equilibration frames, and -kT ln P on the shared grid is then mostly empty bins,
    so a fixed 2 kcal/mol contour has nothing to enclose. A density quantile is defined whatever
    the sample size, which is exactly the property a replicate-agreement check needs.
    """
    from scipy.ndimage import gaussian_filter
    counts, _, _ = np.histogram2d(p[:, 0], p[:, 1], bins=[xe, ye])
    smooth = gaussian_filter(counts, sigma=1.0)
    flat = np.sort(smooth.ravel())[::-1]
    total = flat.sum()
    if total <= 0:
        return None
    level = flat[min(int(np.searchsorted(np.cumsum(flat), fraction * total)), flat.size - 1)]
    if level <= 0:
        return None
    return ax.contour(0.5 * (xe[:-1] + xe[1:]), 0.5 * (ye[:-1] + ye[1:]), smooth.T,
                      levels=[level], colors=[color], linewidths=1.2, linestyles=[ls_])


# --- figure 1: the landscapes -------------------------------------------------------------------
def fig_landscapes(model, proj, xe, ye, fes, fmax, min_count, out: Path) -> Path:
    systems = model["systems"]
    rows, cols = grid_shape(len(systems))
    fig, axes = plt.subplots(rows, cols, figsize=(3.1 * cols + 1.2, 3.0 * rows + 1.0),
                             squeeze=False)
    header(fig, model, "Conformational free-energy landscapes",
           note=f"{len(xe) - 1}x{len(ye) - 1} bins")
    xl, yl = axis_labels(model)
    norm = Normalize(vmin=0, vmax=fmax)
    mesh = None
    for k, ax in enumerate(axes.ravel()):
        if k >= len(systems):
            ax.axis("off")
            continue
        system = systems[k]
        f, counts = fes[system]
        n = int(counts.sum())
        mesh = ax.pcolormesh(xe, ye, np.ma.masked_invalid(f).T, cmap=FES_CMAP, norm=norm,
                             shading="flat", rasterized=True)
        # Contours at whole kcal/mol: the eye reads a basin's shape off lines far better than
        # off a colour gradient, and they survive printing in grey.
        with np.errstate(invalid="ignore"):
            ax.contour(0.5 * (xe[:-1] + xe[1:]), 0.5 * (ye[:-1] + ye[1:]),
                       np.nan_to_num(f.T, nan=np.inf), levels=np.arange(1, fmax + 0.1, 1.0),
                       colors="white", linewidths=0.5, alpha=0.65)
        reps = model["per_system"][system]
        cos = max((max(r["cosine_content"][:2]) for r in reps["replicas"]), default=np.nan)
        flag = "  ⚠ drift" if cos > COSINE_FLAG else ""
        style(ax, xl if k // cols == rows - 1 else "", yl if k % cols == 0 else "",
              f"{system}   {reps['n_replicas']} rep, {n} frames{flag}")
    if mesh is not None:
        cb = fig.colorbar(mesh, ax=axes, fraction=0.02, pad=0.015, extend="max")
        cb.set_label("free energy above the minimum (kcal/mol)", fontsize=7, color=MUTED)
        cb.ax.tick_params(labelsize=6, colors=MUTED)
        cb.outline.set_color(GRID)
    footnote(fig, "Sampling density in energy units, not a converged free energy (D-5): a basin "
                  f"visited twice is a statement about those two visits. Bins with < {min_count} "
                  "frames are left unshaded rather than drawn as a barrier.")
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return out


# --- figure 2: the differences ------------------------------------------------------------------
def fig_differences(model, proj, xe, ye, fes, fmax, min_count, reference: str, out: Path) -> Path | None:
    systems = [s for s in model["systems"] if s != reference]
    if not systems:
        return None
    ref_f, ref_counts = fes[reference]
    rows, cols = grid_shape(len(systems))
    fig, axes = plt.subplots(rows, cols, figsize=(3.1 * cols + 1.2, 3.0 * rows + 1.0),
                             squeeze=False)
    header(fig, model, f"Difference from {reference}")
    xl, yl = axis_labels(model)
    lim = max(1.0, fmax / 2)
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)
    mesh = None
    for k, ax in enumerate(axes.ravel()):
        if k >= len(systems):
            ax.axis("off")
            continue
        system = systems[k]
        f, counts = fes[system]
        # Only where BOTH are sampled. Subtracting a NaN from a number would silently become
        # "this system reaches somewhere the reference does not", which is a claim about the
        # sampling, not about the energetics, and it belongs in the hatching instead.
        both = (counts >= min_count) & (ref_counts >= min_count)
        delta = np.where(both, f - ref_f, np.nan)
        mesh = ax.pcolormesh(xe, ye, np.ma.masked_invalid(delta).T, cmap=DIFF_CMAP, norm=norm,
                             shading="flat", rasterized=True)
        only_here = (counts >= min_count) & (ref_counts < min_count)
        if only_here.any():
            ax.contourf(0.5 * (xe[:-1] + xe[1:]), 0.5 * (ye[:-1] + ye[1:]),
                        only_here.T.astype(float), levels=[0.5, 1.5],
                        colors="none", hatches=["///"])
        style(ax, xl if k // cols == rows - 1 else "", yl if k % cols == 0 else "",
              f"{system} − {reference}")
    if mesh is not None:
        cb = fig.colorbar(mesh, ax=axes, fraction=0.02, pad=0.015, extend="both")
        cb.set_label(f"ΔF (kcal/mol);  blue = deeper than {reference}", fontsize=7, color=MUTED)
        cb.ax.tick_params(labelsize=6, colors=MUTED)
        cb.outline.set_color(GRID)
    footnote(fig, "Shown only where both systems have enough frames to have a free energy; "
                  "hatching marks where this system is sampled and the reference is not, which "
                  "is a difference in sampling and not a measured one.")
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return out


# --- figure 3: whether any of it means anything --------------------------------------------------
def fig_diagnostics(model, proj, xe, ye, out: Path) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.2))
    header(fig, model, "Landscape diagnostics", colorbar=False)
    ax_ts, ax_rep, ax_cos, ax_feat = axes.ravel()

    # (a) implied timescale vs lag -- the justification for the lag that was used
    scan = model.get("timescale_scan_ns")
    if scan:
        lags = np.asarray(scan["lags_ns"], dtype=float)
        ts = np.array([[np.nan if v is None else v for v in row]
                       for row in scan["timescales_ns"]], dtype=float)
        for j in range(ts.shape[1]):
            ax_ts.plot(lags, ts[:, j], marker="o", ms=3.5, lw=1.4,
                       color=REPLICA_COLORS[j % len(REPLICA_COLORS)], label=f"tIC{j + 1}")
        ax_ts.plot(lags, lags, ls=":", lw=1, color=MUTED)
        ax_ts.text(lags[-1], lags[-1], " t = lag", fontsize=6, color=MUTED, va="center")
        ax_ts.axvline(model["lag"]["ns"], color=INK, lw=1, ls="--")
        ax_ts.text(model["lag"]["ns"], ax_ts.get_ylim()[1], " lag used", fontsize=6,
                   color=INK, va="top")
        ax_ts.set_xscale("log")
        ax_ts.set_yscale("log")
        ax_ts.legend(fontsize=6, frameon=False, ncol=2)
        style(ax_ts, "lag (ns)", "implied timescale (ns)",
              "Timescales must stop depending on the lag")
    else:
        empty(ax_ts, "no lag scan — rerun 01_fit_landscape.py --scan")
    ax_ts.set_facecolor(SURFACE)

    # (b) replica-to-replica agreement, as the region holding half of each replica's frames
    fraction = 0.5
    dashes = ["-", "--", ":", "-."]
    coarse_x = np.linspace(xe[0], xe[-1], 29)
    coarse_y = np.linspace(ye[0], ye[-1], 29)
    handles, seen_reps = [], []
    for si, (system, reps) in enumerate(proj.items()):
        for ri, (replica, p) in enumerate(reps):
            hdr_contour(ax_rep, p, coarse_x, coarse_y, fraction,
                        REPLICA_COLORS[ri % len(REPLICA_COLORS)], dashes[si % len(dashes)])
            if replica not in seen_reps:
                seen_reps.append(replica)
                handles.append(plt.Line2D([], [], color=REPLICA_COLORS[ri % len(REPLICA_COLORS)],
                                          lw=1.2, label=replica))
        handles.append(plt.Line2D([], [], color=MUTED, lw=1.2, ls=dashes[si % len(dashes)],
                                  label=system))
    xl, yl = axis_labels(model)
    ax_rep.legend(handles=handles, fontsize=5.5, frameon=False, ncol=2, loc="upper right")
    ax_rep.set_xlim(xe[0], xe[-1])
    ax_rep.set_ylim(ye[0], ye[-1])
    style(ax_rep, xl, yl,
          f"Replicate spread: the region holding {fraction:.0%} of each replica's frames")
    ax_rep.set_facecolor(SURFACE)

    # (c) cosine content -- D-20, applied to the tICs
    rows = [(f"{s}/{r['replica']}", r["cosine_content"][:2])
            for s, v in model["per_system"].items() for r in v["replicas"]]
    if rows:
        y = np.arange(len(rows))
        for j in range(2):
            ax_cos.barh(y + (j - 0.5) * 0.36, [c[j] if len(c) > j else 0 for _, c in rows],
                        height=0.34, color=REPLICA_COLORS[j], label=model["components"][j])
        ax_cos.axvline(COSINE_FLAG, color="#d6452a", lw=1, ls="--")
        ax_cos.text(COSINE_FLAG, -0.9, " not converged →", fontsize=6, color="#d6452a")
        ax_cos.set_yticks(y)
        ax_cos.set_yticklabels([n for n, _ in rows], fontsize=5)
        ax_cos.legend(fontsize=6, frameon=False)
        style(ax_cos, "cosine content", "", "Is the leading mode just drift? (D-20)")
    else:
        empty(ax_cos, "no per-replica record")
    ax_cos.set_facecolor(SURFACE)

    # (d) what the axes are made of
    tica = model.get("tica")
    if tica:
        for j in range(min(2, len(tica["feature_correlations"]))):
            entries = tica["feature_correlations"][j][:8]
            y = np.arange(len(entries)) - (j - 0.5) * 0.36
            ax_feat.barh(y, [e["r"] for e in entries], height=0.34,
                         color=REPLICA_COLORS[j], label=model["components"][j])
            if j == 0:
                ax_feat.set_yticks(np.arange(len(entries)))
                ax_feat.set_yticklabels([e["feature"] for e in entries], fontsize=5.5)
        ax_feat.axvline(0, color=GRID, lw=0.8)
        ax_feat.legend(fontsize=6, frameon=False)
        style(ax_feat, "correlation with the component", "",
              "What the axes are: most-correlated features")
        between = model["between_system_variance_fraction"][:2]
        ax_feat.text(0.99, 0.02,
                     "between-system variance  " + ",  ".join(
                         f"{n} {v:.2f}" for n, v in zip(model["components"], between,
                                                        strict=False)),
                     transform=ax_feat.transAxes, ha="right", fontsize=6, color=MUTED)
    else:
        empty(ax_feat, "named features, not fitted components — nothing to decompose")
    ax_feat.set_facecolor(SURFACE)

    fig.tight_layout(rect=(0, 0.01, 1, 0.94))
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return out


def empty(ax, message: str) -> None:
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=7, style="italic",
            color=MUTED, transform=ax.transAxes, wrap=True)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color(GRID)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--space", choices=("receptor", "ligand", "pair"), default="receptor")
    ap.add_argument("--reference", help="system the difference maps subtract (default: the "
                                        "first ZH853 system, else the first alphabetically)")
    ap.add_argument("--bins", type=int, default=None,
                    help="bins per axis (default: scaled to the sample size)")
    ap.add_argument("--fmax", type=float, default=FMAX, help="colour-scale ceiling, kcal/mol")
    ap.add_argument("--min-count", type=int, default=MIN_COUNT,
                    help="frames per bin below which the bin is left unshaded")
    ap.add_argument("--landscapes", type=Path, default=LAND_ROOT)
    ap.add_argument("--out", type=Path, default=paths.PRODUCT)
    args = ap.parse_args()

    model, proj = load(args.space, args.landscapes)
    xe, ye = edges(proj, auto_bins(proj, args.bins))
    fes = surfaces(proj, xe, ye, args.min_count)
    reference = args.reference or next(
        (s for s in model["systems"] if s.startswith("ZH853")), model["systems"][0])
    if reference not in fes:
        raise SystemExit(f"ERROR: reference {reference} is not among {model['systems']}")

    paths.ensure_dir(args.out)
    stem = f"{PREFIX}_{args.space}"
    stamp = f"{date.today():%Y%m%d}"
    written = [fig_landscapes(model, proj, xe, ye, fes, args.fmax, args.min_count,
                              args.out / f"{stem}_landscapes_{stamp}.png")]
    diff = fig_differences(model, proj, xe, ye, fes, args.fmax, args.min_count, reference,
                           args.out / f"{stem}_differences_{stamp}.png")
    if diff:
        written.append(diff)
    else:
        print("note: only one system — no difference map")
    written.append(fig_diagnostics(model, proj, xe, ye,
                                   args.out / f"{stem}_diagnostics_{stamp}.png"))
    for p in written:
        print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

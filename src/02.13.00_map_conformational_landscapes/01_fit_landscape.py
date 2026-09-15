#!/usr/bin/env python3
"""Fit one shared 2D collective-variable space and project every replica into it.

**Runs locally** (`zh853mor-local`) on the reduced replicas from
[`02.11.00`](../02.11.00_analyze_simulations/README.md). No trajectory, no cluster: the per-frame
feature matrices it needs were written there, because that is the only stage that reads the DCDs.

Three coordinate spaces, chosen with `--space`:

``receptor``  tICA on the pairwise CA-CA distances among the key and functional residues
              (`key_pairs`). Defined for EVERY system, so apo and holo share one landscape.
``ligand``    tICA on every receptor residue's minimum heavy-atom distance to the ligand
              (`min_dist`) -- the Objective-1/2 space. Holo systems only; apo has no ligand.
``pair``      two named features, no fitting: `--x activation:tm3_tm6_ca --y activation:npxxy_tm3_ca`
              is the classic GPCR activation plane. Interpretable by construction, which is the
              thing a tIC has to earn.

**The basis is fitted ONCE, on the pooled data, and every system is projected into it.** Fitting
per system would give each its own axes, and two landscapes drawn on different linear combinations
of distances cannot be laid beside each other however similar they look. `--fit-on <system>`
instead fits on one system and projects the rest, which is the stricter reading -- see below.

**What a pooled fit's "slowest process" actually is.** Systems do not interconvert: no trajectory
turns apo into holo. So a direction that separates two systems has, by construction, an
autocorrelation of ~1 at any lag, and a pooled tICA will happily return it as the slowest mode
with an infinite implied timescale. That is a useful DISCRIMINATIVE axis -- it is the axis along
which the systems differ -- but it is not a relaxation time, and quoting it as one would be wrong.
This step therefore reports, per tIC, the fraction of its variance that lies BETWEEN systems
rather than within them, and re-estimates the timescale separately inside each system. Read the
between-system fraction first: above ~0.5 the tIC is mostly telling you which simulation you are
looking at.

    python 01_fit_landscape.py --space receptor --lag-ns 10 --scan
    python 01_fit_landscape.py --space ligand --systems ZH853_ASP ZH850_ASP ZH831_ASP ZH809_ASP
    python 01_fit_landscape.py --space pair --x activation:tm3_tm6_ca --y activation:npxxy_tm3_ca
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

import numpy as np  # noqa: E402

from zh853mor import convergence as cv  # noqa: E402
from zh853mor import landscape as ls  # noqa: E402
from zh853mor import md, paths  # noqa: E402

ANALYSIS_ROOT = paths.INTERMEDIATE / "02.11.00_analyze_simulations"
OUT_ROOT = paths.INTERMEDIATE / "02.13.00_map_conformational_landscapes"

LAG_NS = 10.0          # tICA lag: long enough to be past the picosecond rattling of a CA distance
N_COMPONENTS = 2       # the landscape is 2D; more are fitted only for the timescale scan
MIN_FRAMES = 200       # below this a replica's post-equilibration window cannot support a fit
BETWEEN_SYSTEM_FLAG = 0.5   # above this, a tIC is a discriminative axis rather than a kinetic one


# --- feature access --------------------------------------------------------------------------
def feature_matrix(result: md.ReplicaResult, space: str) -> tuple[np.ndarray, np.ndarray]:
    """(features, axis) for one replica, post-equilibration. `axis` identifies the columns."""
    arrays = result.arrays()
    if space == "receptor":
        x, axis = arrays.get("key_pairs"), arrays.get("key_pair_ids")
        if x is None or axis is None:
            raise KeyError("no key_pairs in this reduction -- it predates the feature matrices; "
                           "re-run 01_reduce_trajectory.py --force")
    elif space == "ligand":
        x, axis = arrays.get("min_dist"), arrays.get("resids")
    else:
        raise ValueError(f"unknown space {space!r}")
    if x is None or not len(x):
        raise KeyError(f"{space} features are empty for this replica")
    return np.asarray(x, dtype=float)[result.t0:], np.asarray(axis)


def named_feature(result: md.ReplicaResult, name: str) -> np.ndarray:
    """One named per-frame feature, post-equilibration, for `--space pair`.

    Accepts a plain series name from the reduction (`na_d250_dist`, `lig_rgyr`, `thickness`, ...),
    or one of three prefixed forms that reach into the matrices:

        activation:<label>   an activation ruler, by the label 02.11.00 stored
        lig:<resid>          that residue's minimum heavy-atom distance to the ligand
        ca:<a>-<b>           the CA-CA distance between two key residues
    """
    arrays = result.arrays()
    if ":" not in name:
        return result.series(name, equilibrated=True)
    kind, _, rest = name.partition(":")
    if kind == "activation":
        labels = [str(x) for x in arrays["activation_labels"]]
        if rest not in labels:
            raise KeyError(f"activation ruler {rest!r} not among {labels}")
        return np.asarray(arrays["activation"], dtype=float)[result.t0:, labels.index(rest)]
    if kind == "lig":
        resids = [int(r) for r in arrays["resids"]]
        if int(rest) not in resids:
            raise KeyError(f"residue {rest} is not in the construct")
        col = np.asarray(arrays["min_dist"], dtype=float)
        if not len(col):
            raise KeyError("this system is apo: there are no ligand distances")
        return col[result.t0:, resids.index(int(rest))]
    if kind == "ca":
        a, _, b = rest.partition("-")
        pairs = [tuple(int(v) for v in p) for p in arrays["key_pair_ids"]]
        want = (int(a), int(b))
        for candidate in (want, want[::-1]):
            if candidate in pairs:
                return np.asarray(arrays["key_pairs"], dtype=float)[result.t0:,
                                                                    pairs.index(candidate)]
        raise KeyError(f"no CA-CA pair {want} among the key residues")
    raise KeyError(f"unknown feature prefix {kind!r} (expected activation:, lig: or ca:)")


def frames_per_ns(result: md.ReplicaResult) -> float:
    t = result.series("time_ns", equilibrated=False)
    step = float(np.median(np.diff(t))) if t.size > 1 else np.nan
    return 1.0 / step if step and np.isfinite(step) and step > 0 else np.nan


# --- gathering -------------------------------------------------------------------------------
def gather(results: list[md.ReplicaResult], space: str, x: str | None, y: str | None,
           warnings: list[str]) -> tuple[dict[str, list], np.ndarray | None]:
    """{system: [(replica, features, time_ns)]}, dropping what cannot contribute, loudly."""
    out: dict[str, list] = {}
    axis_ref: np.ndarray | None = None
    for r in sorted(results, key=lambda r: (r.system, r.replica)):
        try:
            if space == "pair":
                assert x and y
                columns = [named_feature(r, name) for name in (x, y)]
                # Name the axis that is missing. A feature the reduction recorded as all-NaN --
                # an activation ruler whose CA is absent, the Na+ distance in a build with no
                # D2.50 carboxylate -- otherwise surfaces as "0 usable frames", which reads like
                # a short trajectory rather than a coordinate that does not exist here.
                dead = [name for name, col in zip((x, y), columns, strict=True)
                        if not np.isfinite(col).any()]
                if dead:
                    raise KeyError(f"{', '.join(dead)} is undefined for this system "
                                   "(every frame is NaN)")
                feats = np.column_stack(columns)
                axis = np.array([x, y])
            else:
                feats, axis = feature_matrix(r, space)
        except (KeyError, ValueError) as exc:
            warnings.append(f"{r.system}/{r.replica}: skipped -- {exc}")
            continue
        good = np.isfinite(feats).all(axis=1)
        if good.sum() < MIN_FRAMES:
            warnings.append(f"{r.system}/{r.replica}: only {int(good.sum())} usable "
                            f"post-equilibration frames (< {MIN_FRAMES}); skipped")
            continue
        if good.sum() < feats.shape[0]:
            warnings.append(f"{r.system}/{r.replica}: dropped "
                            f"{feats.shape[0] - int(good.sum())} frames with non-finite features")
        if axis_ref is None:
            axis_ref = axis
        elif axis.shape != axis_ref.shape or not np.array_equal(axis, axis_ref):
            warnings.append(f"{r.system}/{r.replica}: its feature axis differs from the first "
                            "replica's (a different construct or key-residue set); skipped, "
                            "because projecting it would mean something else")
            continue
        time = r.series("time_ns", equilibrated=True)[good]
        out.setdefault(r.system, []).append((r.replica, feats[good], time, r))
    return out, axis_ref


def between_system_fraction(proj: dict[str, list]) -> np.ndarray:
    """Per component, the share of total variance that lies BETWEEN systems, not within them.

    The number that says whether a tIC is a kinetic coordinate or merely a label. Systems never
    interconvert, so a direction separating them has perfect autocorrelation and tICA ranks it
    first; its implied timescale is then an artefact of the comparison, not a property of the
    dynamics.
    """
    groups = [np.vstack([p for _, p, _, _ in reps]) for reps in proj.values()]
    if len(groups) < 2:
        return np.zeros(groups[0].shape[1]) if groups else np.zeros(0)
    allx = np.vstack(groups)
    grand = allx.mean(axis=0)
    total = allx.var(axis=0)
    between = np.average([(g.mean(axis=0) - grand) ** 2 for g in groups],
                         weights=[len(g) for g in groups], axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(total > 0, between / total, 0.0)


def within_system_timescale(series: list[np.ndarray], lag: int, dt_ns: float) -> float:
    """Implied timescale of an already-projected coordinate, inside one system.

    Estimated from the lagged autocorrelation of the projection, pooled over that system's
    replicas -- and never across them, for the same reason the fit does not pair across a replica
    boundary. Unlike the pooled eigenvalue this IS a relaxation time, because within one system
    the coordinate can and does relax.
    """
    num = den = 0.0
    n = 0
    values = [s for s in series if s.size > lag]
    if not values:
        return float("nan")
    grand = np.concatenate(values)
    mean, var = grand.mean(), grand.var()
    if var <= 0:
        return float("nan")
    for s in values:
        a, b = s[:-lag] - mean, s[lag:] - mean
        num += float((a * b).sum())
        n += a.size
    den = var * n
    rho = num / den if den > 0 else np.nan
    if not np.isfinite(rho) or rho <= 0:
        return float("nan")        # decorrelated within the lag: faster than we can resolve
    if rho >= 1:
        return float("inf")        # not decorrelated at all within this window
    return float(-lag * dt_ns / np.log(rho))


def feature_correlations(model: ls.TICA, pooled: np.ndarray, axis: np.ndarray,
                         top: int = 8) -> list[list[dict]]:
    """The `top` features most correlated with each tIC -- what the axis actually MEANS.

    Correlation rather than the raw loading: a loading is in whitened units and a feature with a
    large loading and almost no variance contributes nothing you could point at in a structure.
    """
    proj = model.transform(pooled)
    out = []
    for k in range(proj.shape[1]):
        p = proj[:, k]
        pc = p - p.mean()
        f = pooled - pooled.mean(axis=0)
        denom = np.sqrt((pc ** 2).sum() * (f ** 2).sum(axis=0))
        with np.errstate(divide="ignore", invalid="ignore"):
            corr = np.where(denom > 0, (f * pc[:, None]).sum(axis=0) / denom, 0.0)
        order = np.argsort(np.abs(corr))[::-1][:top]
        out.append([{"feature": _axis_label(axis, int(i)), "r": round(float(corr[i]), 3)}
                    for i in order])
    return out


def _axis_label(axis: np.ndarray, i: int) -> str:
    entry = axis[i]
    if np.ndim(entry) and len(entry) == 2:
        return f"CA {int(entry[0])}-{int(entry[1])}"
    try:
        return f"residue {int(entry)}"
    except (TypeError, ValueError):
        return str(entry)


# --- the fit ----------------------------------------------------------------------------------
def run(args) -> dict:
    results = md.load_replicas(args.analysis)
    if not results:
        raise SystemExit(f"ERROR: no reduced replicas under {args.analysis}. Run "
                         "02.11.00's 01_reduce_trajectory.py on the cluster and copy them back.")
    if args.systems:
        results = [r for r in results if r.system in set(args.systems)]
        if not results:
            raise SystemExit(f"ERROR: none of {args.systems} are among the reduced systems: "
                             f"{sorted({r.system for r in md.load_replicas(args.analysis)})}")

    warnings: list[str] = []
    proj_in, axis = gather(results, args.space, args.x, args.y, warnings)
    if not proj_in:
        raise SystemExit("ERROR: no replica could contribute features.\n  "
                         + "\n  ".join(warnings))
    if len(proj_in) < 2:
        warnings.append(f"only one system ({next(iter(proj_in))}) has usable features; this is a "
                        "landscape, not a comparison")

    rates = {fpn for reps in proj_in.values() for *_, r in reps
             if np.isfinite(fpn := frames_per_ns(r))}
    if len({round(v, 6) for v in rates}) > 1:
        warnings.append(f"replicas were written at different frame rates ({sorted(rates)} per ns); "
                        "the lag is set from the fastest, so it is not the same simulated time "
                        "for every replica")
    fpn = max(rates) if rates else 10.0
    lag = max(1, int(round(args.lag_ns * fpn)))
    dt_ns = 1.0 / fpn

    model = None
    scan = None
    if args.space == "pair":
        projected = {s: [(name, f, t, r) for name, f, t, r in reps] for s, reps in proj_in.items()}
        components = [args.x, args.y]
    else:
        if args.fit_on:
            if args.fit_on not in proj_in:
                raise SystemExit(f"ERROR: --fit-on {args.fit_on} is not among "
                                 f"{sorted(proj_in)}")
            fit_set = [f for _, f, _, _ in proj_in[args.fit_on]]
            warnings.append(f"basis fitted on {args.fit_on} alone and applied to the rest: the "
                            "axes are that system's kinetic modes, and another system's density "
                            "is being read in coordinates chosen without it")
        else:
            fit_set = [f for reps in proj_in.values() for _, f, _, _ in reps]
        if args.transform == "inverse":
            # Contact-like: 1/d compresses the long distances, where a 2 A change means nothing,
            # and expands the short ones, where it is a contact forming or breaking.
            fit_set = [1.0 / np.clip(f, 1e-3, None) for f in fit_set]
        model = ls.fit_tica(fit_set, lag=lag, n_components=args.components,
                            kinetic_map=not args.no_kinetic_map)
        if args.scan:
            lags = sorted({max(1, int(round(f * fpn)))
                           for f in (0.5, 1, 2, 5, 10, 20, 50, 100) if f <= args.lag_ns * 10})
            scan = ls.implied_timescale_scan(
                fit_set, lags=lags, n_components=min(4, args.components + 2),
                kinetic_map=not args.no_kinetic_map)
        projected = {}
        for system, reps in proj_in.items():
            rows = []
            for name, f, t, r in reps:
                feats = 1.0 / np.clip(f, 1e-3, None) if args.transform == "inverse" else f
                rows.append((name, model.transform(feats), t, r))
            projected[system] = rows
        components = [f"tIC{i + 1}" for i in range(args.components)]

    # --- what the axes are worth ---------------------------------------------------------------
    between = between_system_fraction(projected)
    per_system = {}
    for system, reps in projected.items():
        per_system[system] = {
            "n_replicas": len(reps),
            "n_frames": int(sum(len(p) for _, p, _, _ in reps)),
            "replicas": [{
                "replica": name,
                "n_frames": int(len(p)),
                "t0_ns": round(float(r.t0_ns), 2),
                # The D-20 guard, in tICA form: on a coordinate that is still drifting, the
                # projection is the half-cosine of free diffusion and the landscape is a picture
                # of the drift rather than of an ensemble.
                "cosine_content": [round(float(cv.cosine_content(p[:, k], index=1)), 3)
                                   for k in range(p.shape[1])],
            } for name, p, _, r in reps],
            "timescale_ns": [round(within_system_timescale([p[:, k] for _, p, _, _ in reps],
                                                           lag, dt_ns), 2)
                             for k in range(len(components))],
        }

    report: dict = {
        "space": args.space,
        "components": components,
        "systems": sorted(projected),
        "lag": {"frames": lag, "ns": round(lag * dt_ns, 3), "requested_ns": args.lag_ns},
        "fit_on": args.fit_on or ("pooled" if args.space != "pair"
                                  else "no fit — named features"),
        "between_system_variance_fraction": [round(float(v), 3) for v in between],
        "between_system_flag": BETWEEN_SYSTEM_FLAG,
        "per_system": per_system,
        "warnings": warnings,
    }
    if model is not None:
        pooled = np.vstack([f for reps in proj_in.values() for _, f, _, _ in reps])
        if args.transform == "inverse":
            pooled = 1.0 / np.clip(pooled, 1e-3, None)
        report["tica"] = {
            "eigenvalues": [round(float(v), 4) for v in model.eigenvalues],
            "pooled_timescale_ns": [None if not np.isfinite(v) else round(float(v), 2)
                                    for v in model.timescales(dt_ns)],
            "rank": model.rank, "n_pairs": model.n_samples,
            "n_features": int(pooled.shape[1]),
            "kinetic_map": not args.no_kinetic_map, "transform": args.transform,
            "feature_correlations": feature_correlations(model, pooled, axis),
        }

    out_dir = paths.ensure_dir(args.out / args.space)
    payload = {}
    for system, reps in projected.items():
        for name, p, t, _ in reps:
            payload[f"proj__{system}__{name}"] = p.astype(np.float32)
            payload[f"time__{system}__{name}"] = t.astype(np.float32)
    if scan is not None:
        payload["scan_lags"], payload["scan_timescales"] = scan[0], scan[1]
        report["timescale_scan_ns"] = {
            "lags_ns": [round(float(v * dt_ns), 3) for v in scan[0]],
            "timescales_ns": [[None if not np.isfinite(v) else round(float(v * dt_ns), 2)
                               for v in row] for row in scan[1]],
        }
    np.savez_compressed(out_dir / "projections.npz", **payload)
    (out_dir / "model.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def report_lines(report: dict) -> list[str]:
    """The terminal report and the product markdown are the SAME text, built once.

    They are the record of what basis a landscape was drawn on, and two formattings of it would
    eventually disagree about a number that matters (which lag, fitted on what, how much of tIC1
    is between systems).
    """
    out: list[str] = []
    add = out.append
    add(f"\nspace: {report['space']}   basis: {report['fit_on']}   "
          f"lag: {report['lag']['ns']} ns ({report['lag']['frames']} frames)")
    tica = report.get("tica")
    if tica:
        add(f"  {tica['n_features']} features -> rank {tica['rank']} after whitening, "
              f"{tica['n_pairs']} lagged pairs")
        for k, name in enumerate(report["components"]):
            ts = tica["pooled_timescale_ns"][k]
            add(f"  {name}: lambda={tica['eigenvalues'][k]:.3f}  "
                  f"pooled timescale={'>window' if ts is None else f'{ts:g} ns'}  "
                  f"between-system variance={report['between_system_variance_fraction'][k]:.2f}")
            tops = ", ".join(f"{c['feature']} ({c['r']:+.2f})"
                             for c in tica["feature_correlations"][k][:5])
            add(f"        most correlated features: {tops}")
    flagged = [n for n, v in zip(report["components"],
                                 report["between_system_variance_fraction"], strict=False)
               if v > report["between_system_flag"]]
    if flagged:
        add(f"  NOTE: {', '.join(flagged)} carry more variance BETWEEN systems than within. "
              "They separate the simulations, which is what makes the comparison legible, but "
              "their pooled timescales are not relaxation times.")

    add(f"\n  {'system':<14}{'reps':>5}{'frames':>8}   " +
          "  ".join(f"{c} ts/ns" for c in report["components"]) + "   max cosine")
    for system, s in report["per_system"].items():
        cos = max((max(r["cosine_content"]) for r in s["replicas"]), default=float("nan"))
        ts = "  ".join(f"{v:>9.1f}" if np.isfinite(v) else f"{'--':>9}"
                       for v in s["timescale_ns"])
        mark = "  <-- drift?" if cos > 0.5 else ""
        add(f"  {system:<14}{s['n_replicas']:>5}{s['n_frames']:>8}   {ts}   {cos:>8.2f}{mark}")
    for w in report["warnings"]:
        add(f"  WARNING: {w}")
    return out


def write_markdown(report: dict, path: Path) -> Path:
    body = "\n".join(report_lines(report))
    path.write_text(
        f"# 02.13.00 — conformational landscape ({report['space']} space)\n\n"
        f"Generated {date.today():%Y-%m-%d} by `01_fit_landscape.py`. Figures: "
        f"`02.13.00_{report['space']}_*.png`.\n\n"
        "```\n" + body.strip("\n") + "\n```\n\n"
        "**How to read this.** The landscape is a sampling density in energy units, not a "
        "converged free energy (D-5). Two systems are comparable only because the basis was "
        "fitted once and every system projected into it; a component whose between-system "
        "variance fraction exceeds "
        f"{report['between_system_flag']} separates the simulations rather than describing a "
        "relaxation, so its pooled timescale is not a rate. A replica whose cosine content "
        "exceeds 0.5 is still drifting, and its density is a picture of that drift (D-20).\n")
    return path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--space", choices=("receptor", "ligand", "pair"), default="receptor")
    ap.add_argument("--x", help="[--space pair] first feature, e.g. activation:tm3_tm6_ca")
    ap.add_argument("--y", help="[--space pair] second feature, e.g. lig:299")
    ap.add_argument("--systems", nargs="+", help="restrict the comparison set")
    ap.add_argument("--fit-on", help="fit the basis on this system alone, then project the rest")
    ap.add_argument("--lag-ns", type=float, default=LAG_NS, help=f"tICA lag (default {LAG_NS})")
    ap.add_argument("--components", type=int, default=N_COMPONENTS)
    ap.add_argument("--transform", choices=("none", "inverse"), default="none",
                    help="'inverse' uses 1/d, which weights contact formation over bulk distance")
    ap.add_argument("--no-kinetic-map", action="store_true",
                    help="do not scale each tIC by its eigenvalue")
    ap.add_argument("--scan", action="store_true",
                    help="also scan implied timescales against lag (the lag justification)")
    ap.add_argument("--analysis", type=Path, default=ANALYSIS_ROOT)
    ap.add_argument("--out", type=Path, default=OUT_ROOT)
    ap.add_argument("--product", type=Path, default=paths.PRODUCT,
                    help="where the markdown record of the fit is written")
    args = ap.parse_args()

    if args.space == "pair" and not (args.x and args.y):
        raise SystemExit("ERROR: --space pair needs both --x and --y, e.g. "
                         "--x activation:tm3_tm6_ca --y activation:npxxy_tm3_ca")
    if args.space != "pair" and (args.x or args.y):
        raise SystemExit(f"ERROR: --x/--y only apply to --space pair, not {args.space}")

    report = run(args)
    for line in report_lines(report):
        print(line)
    md = write_markdown(report, paths.ensure_dir(args.product) /
                        f"02.13.00_{args.space}_model_{date.today():%Y%m%d}.md")
    print(f"\n  -> {args.out / args.space}/  (projections.npz, model.json)")
    print(f"  -> {md}")
    print("  then: python 02_figures.py --space " + args.space)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

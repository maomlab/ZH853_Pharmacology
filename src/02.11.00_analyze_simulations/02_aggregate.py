#!/usr/bin/env python3
"""Aggregate the reduced replicas into QC verdicts, tables and a report.

**Runs locally** (`zh853mor-local`) on what `01_reduce_trajectory.py` wrote -- a few hundred kB
per replica -- so no trajectory has to leave the cluster.

Two questions, in this order:

1. **Is the sampling good enough to believe?** Thermodynamics, membrane, receptor and ligand
   stability, and -- the part a plateau cannot answer -- effective sample size after
   equilibration, residual drift, agreement between replicas (R-hat) and the cosine content of
   the leading principal components. Verdicts are PASS/WARN/FAIL against the thresholds below.
2. **What do the simulations say?** Contact occupancy with replicate error bars for the anchors
   and the candidate ZH853-distinctive contacts (Objectives 1-2), direct vs water-mediated
   H6.52, ligand pose retention, the activation rulers with and without ligand, and the D2.50
   ASP/ASH comparison (SPECIFICATION D-11).

Run: ``python src/02.11.00_analyze_simulations/02_aggregate.py``  (or ``make sim-aggregate``).
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from datetime import date
from pathlib import Path

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

import numpy as np  # noqa: E402

from zh853mor import convergence as cv  # noqa: E402
from zh853mor import md, paths, structure  # noqa: E402

ANALYSIS_ROOT = paths.INTERMEDIATE / "02.11.00_analyze_simulations"
PREFIX = "02.11.00"

# --- QC thresholds ----------------------------------------------------------------------------
# Physical ranges follow `check_equilibration.py` (9:1 POPC:cholesterol, OPC, 0.15 M NaCl, 310 K).
# The structural tolerances are deliberately LOOSER here: that script judges a 2.25 ns restrained
# ramp, where the receptor should not have moved, while these are 500 ns of free dynamics, where
# it should. A production Ca RMSD of 2-3 A against the starting model is a relaxing protein, not
# a broken one; 3.5 A is where it stops being a refinement and starts being a different structure.
TEMP_TARGET_K, TEMP_TOL_K = 310.0, 2.0
DENSITY_RANGE = (0.98, 1.10)          # g/mL
APL_NET_RANGE = (55.0, 70.0)          # A^2/lipid, protein cross-section removed (lower bound)
THICKNESS_RANGE = (34.0, 44.0)        # A, phosphate-to-phosphate
CA_RMSD_MAX = 3.5                     # A, whole receptor vs the staged OPM-oriented receptor
TM_RMSD_MAX = 2.5                     # A, membrane-embedded Ca only
LIG_RMSD_MAX = 4.0                    # A, receptor-aligned ligand vs the deposited pose
SS_RANGE = (1.90, 2.30)               # A, C142-C219
REGISTRATION_DRIFT_A = 2.0            # A of vertical drift out of the OPM slab
DENSITY_DRIFT_PER_NS = 1e-4           # g/mL/ns; 500 ns x this is 0.05, well inside DENSITY_RANGE
APL_DRIFT_PER_NS = 0.02               # A^2/lipid/ns: a bilayer still condensing is not sampled yet

# Sampling-quality thresholds.
MIN_EFFECTIVE_SAMPLES = 20            # per replica, after equilibration; below this a mean is
                                      # a single observation with an error bar drawn around it
RHAT_MAX = 1.2                        # replica agreement; 1.1 is conventional, relaxed for n=3
COSINE_MAX = 0.5                      # Hess: above this, PC1 is diffusion, not sampling
NA_SITE_CUT = 3.2                     # A, Na+ to a D2.50 carboxylate O = occupying the site
DRIFT_SIGMA = 2.0                     # |slope| > this many standard errors = a resolved drift

# Observables carried through the convergence analysis. Chosen because each is load-bearing for a
# claim: the first three for structural stability, the last two for the membrane environment.
CONVERGENCE_OBSERVABLES = ["rmsd_ca_tm", "rmsd_ca_all", "lig_rmsd_pose", "apl_net", "thickness"]


@dataclass
class Check:
    name: str
    value: str
    verdict: str
    detail: str


def judge(name: str, value: float, lo: float, hi: float, fmt: str = "{:.2f}",
          unit: str = "", soft: bool = False, detail: str = "") -> Check:
    if not np.isfinite(value):
        return Check(name, "n/a", "SKIP", detail or "not measured")
    ok = lo <= value <= hi
    return Check(name, fmt.format(value) + unit, "PASS" if ok else ("WARN" if soft else "FAIL"),
                 detail or f"expected {fmt.format(lo)}-{fmt.format(hi)}{unit}")


def combine(values: list[float]) -> tuple[float, float]:
    """(mean, standard error) across replicas -- the replicate spread the QC section asks for."""
    arr = np.array([v for v in values if np.isfinite(v)], dtype=float)
    if arr.size == 0:
        return float("nan"), float("nan")
    if arr.size == 1:
        return float(arr[0]), float("nan")
    return float(arr.mean()), float(arr.std(ddof=1) / np.sqrt(arr.size))


def pm(mean: float, sem: float, digits: int = 2) -> str:
    """"1.83 +- 0.07", or "1.83 (1 rep)" when there is no spread to quote."""
    if not np.isfinite(mean):
        return "n/a"
    if not np.isfinite(sem):
        return f"{mean:.{digits}f} (1 rep)"
    return f"{mean:.{digits}f} +- {sem:.{digits}f}"


def table(headers: list[str], rows: list[list[str]]) -> list[str]:
    """A GitHub-flavoured markdown table. Cell pipes are escaped: several thresholds are written
    as |drift| <= x, which would otherwise split the cell and shift every column after it."""
    def cells(row: list[str]) -> str:
        return "| " + " | ".join(c.replace("|", "\\|") for c in row) + " |"

    return [cells(headers), "|" + "|".join("---" for _ in headers) + "|",
            *[cells(r) for r in rows]]


def system_qc(reps: list[md.ReplicaResult]) -> list[Check]:
    """Every production QC check for one system, averaged over its replicas."""
    def mean_of(name: str) -> tuple[float, float]:
        return combine([r.mean_of(name) for r in reps])

    def drift_of(name: str) -> tuple[float, float]:
        slopes = [r.summary.get("means", {}).get(name, {}).get("drift_per_ns", np.nan)
                  for r in reps]
        return combine([float(s) for s in slopes])

    checks: list[Check] = []
    t, _ = mean_of("log_temperature")
    checks.append(judge("temperature", t, TEMP_TARGET_K - TEMP_TOL_K, TEMP_TARGET_K + TEMP_TOL_K,
                        unit=" K", detail=f"{TEMP_TARGET_K:.0f} +- {TEMP_TOL_K:.0f} K"))
    d, _ = mean_of("log_density")
    checks.append(judge("density", d, *DENSITY_RANGE, fmt="{:.3f}", unit=" g/mL"))
    dd, _ = drift_of("log_density")
    checks.append(judge("density drift", abs(dd), 0.0, DENSITY_DRIFT_PER_NS, fmt="{:.2e}",
                        unit=" g/mL/ns", soft=True,
                        detail=f"|drift| <= {DENSITY_DRIFT_PER_NS:.0e} g/mL/ns"))
    a, _ = mean_of("apl_net")
    checks.append(judge("area per lipid (net)", a, *APL_NET_RANGE, unit=" A^2",
                        soft=True, detail=f"{APL_NET_RANGE[0]:.0f}-{APL_NET_RANGE[1]:.0f} A^2; "
                                          "hull-corrected, so a lower bound"))
    ad, _ = drift_of("apl_net")
    checks.append(judge("area per lipid drift", abs(ad), 0.0, APL_DRIFT_PER_NS, fmt="{:.3f}",
                        unit=" A^2/ns", soft=True,
                        detail=f"|drift| <= {APL_DRIFT_PER_NS} A^2/lipid/ns"))
    th, _ = mean_of("thickness")
    checks.append(judge("bilayer thickness", th, *THICKNESS_RANGE, unit=" A"))
    reg, _ = mean_of("registration_z")
    checks.append(judge("OPM registration", abs(reg), 0.0, REGISTRATION_DRIFT_A, unit=" A",
                        soft=True, detail=f"|TM centre - midplane| <= {REGISTRATION_DRIFT_A} A"))
    rall, _ = mean_of("rmsd_ca_all")
    checks.append(judge("Ca RMSD (whole)", rall, 0.0, CA_RMSD_MAX, unit=" A",
                        detail=f"<= {CA_RMSD_MAX} A vs the staged receptor"))
    rtm, _ = mean_of("rmsd_ca_tm")
    checks.append(judge("Ca RMSD (TM)", rtm, 0.0, TM_RMSD_MAX, unit=" A",
                        detail=f"<= {TM_RMSD_MAX} A"))
    ss, _ = mean_of("ss_dist")
    checks.append(judge("C142-C219 disulfide", ss, *SS_RANGE, unit=" A"))
    if any(r.summary.get("n_ligand_atoms") for r in reps):
        lig, _ = mean_of("lig_rmsd_pose")
        checks.append(judge("ligand RMSD vs pose", lig, 0.0, LIG_RMSD_MAX, unit=" A",
                            detail=f"<= {LIG_RMSD_MAX} A, receptor-aligned"))
    n_eff = min(float(r.summary.get("equilibration", {}).get("n_eff", np.nan)) for r in reps)
    checks.append(judge("effective samples (worst replica)", n_eff, MIN_EFFECTIVE_SAMPLES, 1e9,
                        fmt="{:.0f}", detail=f">= {MIN_EFFECTIVE_SAMPLES} after equilibration"))
    return checks


def convergence_rows(system: str, reps: list[md.ReplicaResult]) -> list[list[str]]:
    """R-hat and blocking for each convergence observable, plus PC1 cosine content."""
    rows = []
    for name in CONVERGENCE_OBSERVABLES:
        chains = []
        for r in reps:
            try:
                s = r.series(name)
            except KeyError:
                continue
            s = s[np.isfinite(s)]
            if s.size > 10:
                chains.append(s)
        if not chains:
            continue
        rhat = cv.gelman_rubin(chains) if len(chains) > 1 else float("nan")
        gs = [cv.statistical_inefficiency(c) for c in chains]
        block = []
        for c in chains:
            _, sems = cv.block_sem_curve(c)
            # The LAST few points are the plateau, i.e. the error bar that accounts for the
            # correlation; the first point is the naive std/sqrt(N) that ignores it.
            block.append(sems[-3:].mean() if sems.size >= 3 else np.nan)
        rows.append([system, name, f"{np.mean(gs):.1f}",
                     f"{np.mean([c.size for c in chains]) / np.mean(gs):.0f}",
                     "n/a" if not np.isfinite(rhat) else f"{rhat:.2f}",
                     f"{np.nanmean(block):.3f}",
                     "" if not np.isfinite(rhat) else ("ok" if rhat <= RHAT_MAX else "REPLICAS DISAGREE")])
    cosines = [r.summary.get("pca", {}).get("cosine_content", [np.nan])[0] for r in reps]
    c_mean, c_sem = combine([float(c) for c in cosines])
    rows.append([system, "PC1 cosine content", "-", "-", "-", pm(c_mean, c_sem),
                 "ok" if c_mean <= COSINE_MAX else "PC1 IS DIFFUSION"])
    return rows


def occupancy_table(systems: dict[str, list[md.ReplicaResult]]) -> tuple[list[str], list[list[str]]]:
    """Anchor contact occupancy (fraction of frames within 4.5 A), mean +- replicate SEM."""
    names = [s for s in systems if any(r.summary.get("n_ligand_atoms") for r in systems[s])]
    headers = ["residue", "BW", *names]
    rows = []
    for resid, label in md.ANCHORS.items():
        row = [f"{resid}", label]
        for s in names:
            vals = [float(r.summary.get("occupancy", {}).get(str(resid), np.nan))
                    for r in systems[s]]
            row.append(pm(*combine(vals)))
        rows.append(row)
    return headers, rows


def water_bridge_table(systems: dict[str, list[md.ReplicaResult]]) -> list[list[str]]:
    names = [s for s in systems if any(r.summary.get("n_ligand_atoms") for r in systems[s])]
    rows = []
    for resid, label in md.ANCHORS.items():
        row = [f"{resid}", label]
        for s in names:
            direct = combine([float(r.summary.get("polar_occupancy", {}).get(str(resid), np.nan))
                              for r in systems[s]])
            bridged = combine([float(r.summary.get("water_bridges", {}).get(str(resid), np.nan))
                               for r in systems[s]])
            row.append(f"{pm(*direct)} / {pm(*bridged)}")
        rows.append(row)
    return rows


def distinctive_contacts(systems: dict[str, list[md.ReplicaResult]], top: int = 12
                         ) -> list[list[str]]:
    """Residues whose occupancy in the ZH853 system most exceeds every other ligand's.

    The MD counterpart of the static Objective-1 analysis: a contact is ZH853-distinctive only if
    it is both persistent in ZH853 and absent in the analogs, which a single structure cannot say.
    """
    per_system: dict[str, np.ndarray] = {}
    resids: np.ndarray | None = None
    for s, reps in systems.items():
        if not any(r.summary.get("n_ligand_atoms") for r in reps):
            continue
        stacks = []
        for r in reps:
            arrays = r.arrays()
            if "occupancy" not in arrays or not arrays["occupancy"].size:
                continue
            stacks.append(arrays["occupancy"])
            resids = arrays["resids"]
        if stacks:
            per_system[s] = np.mean(stacks, axis=0)
    zh = [s for s in per_system if s.startswith("ZH853")]
    others = [s for s in per_system if not s.startswith("ZH853")]
    if not zh or not others or resids is None:
        return []
    ref = np.mean([per_system[s] for s in zh], axis=0)
    rival = np.max([per_system[s] for s in others], axis=0)
    delta = ref - rival
    order = np.argsort(delta)[::-1][:top]
    return [[f"{int(resids[i])}", structure.bw(int(resids[i])), f"{ref[i]:.2f}",
             f"{rival[i]:.2f}", f"{delta[i]:+.2f}"] for i in order if delta[i] > 0]


def activation_table(systems: dict[str, list[md.ReplicaResult]]) -> list[list[str]]:
    rows = []
    for s, reps in systems.items():
        row = [s, f"{len(reps)}"]
        for label, _, _ in md.ACTIVATION_DISTANCES:
            row.append(pm(*combine([r.mean_of(label) for r in reps])))
        na_occ = []
        for r in reps:
            try:
                d = r.series("na_d250_dist")
            except KeyError:
                continue
            d = d[np.isfinite(d)]
            if d.size:
                na_occ.append(float((d <= NA_SITE_CUT).mean()))
        row.append(pm(*combine(na_occ)))
        rows.append(row)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--analysis-root", type=Path, default=ANALYSIS_ROOT)
    ap.add_argument("--date", default=f"{date.today():%Y%m%d}", help="output date code")
    args = ap.parse_args()

    results = md.load_replicas(args.analysis_root)
    if not results:
        raise SystemExit(
            f"ERROR: no reduced replicas under {args.analysis_root}.\n"
            "  Run 01_reduce_trajectory.py on the cluster first, then copy\n"
            "  intermediate/02.11.00_analyze_simulations/ back (see this stage's README).")
    systems = md.group_by_system(results)
    paths.ensure_dir(paths.PRODUCT)

    # --- tables ------------------------------------------------------------------------------
    qc_rows, all_checks = [], {}
    for s, reps in systems.items():
        checks = system_qc(reps)
        all_checks[s] = checks
        for c in checks:
            qc_rows.append([s, c.name, c.value, c.verdict, c.detail])
    qc_csv = paths.PRODUCT / f"{PREFIX}_simulation_qc_{args.date}.csv"
    qc_csv.write_text("system,check,value,verdict,expected\n" + "".join(
        ",".join(f'"{f}"' for f in row) + "\n" for row in qc_rows))

    occ_headers, occ_rows = occupancy_table(systems)
    occ_csv = paths.PRODUCT / f"{PREFIX}_contact_occupancy_{args.date}.csv"
    occ_csv.write_text(",".join(occ_headers) + "\n" + "".join(
        ",".join(f'"{f}"' for f in row) + "\n" for row in occ_rows))

    # --- report ------------------------------------------------------------------------------
    total_ns = sum(float(r.summary.get("length_ns", 0.0)) for r in results)
    lines = [
        f"# Production MD: quality control and interaction analysis ({args.date})",
        "",
        f"Generated by `src/02.11.00_analyze_simulations/02_aggregate.py` from "
        f"{len(results)} replicas of {len(systems)} system(s), {total_ns / 1000:.2f} us in total.",
        "Residue numbers are human OPRM1 (P35372); occupancy is the fraction of post-equilibration",
        f"frames with a heavy-atom contact within {md.CONTACT_CUT} A, the same criterion as the",
        "static cryo-EM fingerprint in `03.01.00`, so the two are directly comparable.",
        "",
        "## 1. Inventory and sampling",
        "",
    ]
    inv = []
    for s, reps in systems.items():
        lengths = [float(r.summary.get("length_ns", np.nan)) for r in reps]
        t0s = [float(r.summary.get("equilibration", {}).get("t0_ns", np.nan)) for r in reps]
        neff = [float(r.summary.get("equilibration", {}).get("n_eff", np.nan)) for r in reps]
        # Against the length the build was configured for: "143" and "143 of 500" are very
        # different statements about a replica, and only the second is checkable.
        targets = [float(r.summary.get("target_ns") or np.nan) for r in reps]
        target = np.nanmean(targets) if np.isfinite(targets).any() else np.nan
        per_replica = f"{np.nanmean(lengths):.0f}"
        if np.isfinite(target) and target > 0:
            per_replica += f" / {target:.0f} ({100 * np.nanmean(lengths) / target:.0f}%)"
        inv.append([s, f"{len(reps)}", f"{np.nansum(lengths):.0f}",
                    per_replica, f"{np.nanmean(t0s):.0f}",
                    f"{np.nanmin(neff):.0f}",
                    reps[0].summary.get("ligand_reference", "-")])
    lines += table(["system", "replicas", "total ns", "ns/replica (of target)", "mean t0 (ns)",
                    "min n_eff", "ligand reference"], inv)
    lines += [
        "",
        "`t0` is where `convergence.detect_equilibration` places the end of the relaxation on the",
        "TM Ca RMSD; every mean below is taken after it. `n_eff` is the number of INDEPENDENT",
        "samples that leaves -- not the number of frames.",
        "",
    ]
    # Anything the reduction noticed about a replica belongs in the report, not only in its JSON:
    # a ligand whose pose reference could not be resolved, or a tautomer that disagrees with the
    # prepared receptor, changes how the numbers below should be read.
    flagged: dict[str, list[str]] = {}
    for r in results:
        for w in r.summary.get("warnings", []):
            flagged.setdefault(w, []).append(f"{r.system}/{r.replica}")
    if flagged:
        lines += ["### Notes from the reduction", ""]
        # Grouped by note: the same condition usually holds for every replica of a system, and a
        # row per replica would bury the one note that differs.
        lines += table(["note", "where"], [
            [note, f"{len(who)} replica(s): " + ", ".join(who[:4])
             + (f", +{len(who) - 4} more" if len(who) > 4 else "")]
            for note, who in flagged.items()])
        lines += [""]
    lines += ["## 2. Quality control", ""]
    lines += table(["system", "check", "value", "verdict", "expected"], qc_rows)
    failed = {s: [c.name for c in cs if c.verdict == "FAIL"] for s, cs in all_checks.items()}
    warned = {s: [c.name for c in cs if c.verdict == "WARN"] for s, cs in all_checks.items()}
    lines += ["", "**Verdict.** " + (
        "; ".join(f"{s}: FAIL ({', '.join(v)})" for s, v in failed.items() if v)
        or "no FAIL in any system") + ".",
        "Warnings: " + ("; ".join(f"{s}: {', '.join(v)}" for s, v in warned.items() if v)
                        or "none") + ".",
        "",
        "## 3. Convergence beyond the plateau",
        "",
    ]
    conv_rows = [row for s, reps in systems.items() for row in convergence_rows(s, reps)]
    lines += table(["system", "observable", "g (frames/sample)", "n_eff", "R-hat",
                    "blocked SEM", "flag"], conv_rows)
    lines += [
        "",
        f"R-hat > {RHAT_MAX} means the replicas are not sampling the same distribution -- with",
        "three replicas this is a flag, not a measurement. A PC1 cosine content above",
        f"{COSINE_MAX} means the leading 'collective motion' is indistinguishable from free",
        "diffusion (Hess 2002), i.e. the run is too short for the motion it appears to show.",
        "",
        "## 4. Interactions (Objectives 1-2)",
        "",
    ]
    if occ_rows:
        lines += table(occ_headers, occ_rows)
        lines += [
            "",
            "Occupancy of the anchor contacts, mean +- SEM over replicas. Direct polar contact",
            f"(N/O-N/O within {md.HBOND_CUT} A) vs bridging waters (count within "
            f"{md.HBOND_CUT} A of both partners):",
            "",
        ]
        lines += table(occ_headers, water_bridge_table(systems))
        lines += [
            "",
            "H6.52 (His299) is the position to read carefully: the deposited pose puts it at",
            "4.89 A, outside the direct shell, and the canonical MOR contact at this position is",
            "water-mediated. A low direct occupancy with a non-zero bridge count is that",
            "interaction, not its absence.",
            "",
        ]
        ligand_systems = [s for s in systems
                          if any(r.summary.get("n_ligand_atoms") for r in systems[s])]
        comparable = (any(s.startswith("ZH853") for s in ligand_systems)
                      and any(not s.startswith("ZH853") for s in ligand_systems))
        dist_rows = distinctive_contacts(systems)
        if dist_rows:
            lines += ["### ZH853-distinctive contacts", "",
                      "Occupancy in the ZH853 system minus the highest occupancy among the other",
                      "ligand systems. Positive = held by ZH853 and not by the analogs, i.e. the",
                      "mutation candidates that should spare the comparators (Objective 2).", ""]
            lines += table(["residue", "BW", "ZH853", "best analog", "delta"], dist_rows)
            lines += [""]
        elif comparable:
            lines += ["*No residue is more occupied in ZH853 than in every analog: on this",
                      "sampling the analogs hold the same contacts, which is itself the Objective-2",
                      "answer -- no contact here discriminates ZH853 from them.*", ""]
        else:
            lines += ["*ZH853-distinctive contacts need the ZH853 system AND at least one analog",
                      "system; only " + ", ".join(ligand_systems) + " reduced so far.*", ""]
    else:
        lines += ["*No ligand-bearing system has been reduced yet: every system present is apo,*",
                  "*so contact occupancy is undefined.*", ""]

    lines += ["## 5. Activation state and the D2.50 sodium site", "",
              *table(["system", "replicas", "R3.50-T6.34 (A)", "R3.50-Y7.53 (A)",
                      f"Na+ at D2.50 (<{NA_SITE_CUT} A)"], activation_table(systems)),
              "",
              "The systems carry no transducer (SPECIFICATION D-10 system B), so the TM6 rulers",
              "measure whether the ligand alone holds the active state: a contraction of",
              "R3.50-T6.34 over the run is the receptor relaxing towards inactive, and comparing",
              "apo against the holo systems is the test of whether the ligand prevents it.",
              "",
              "The sodium column is the direct test of D-11: D2.50 was built both charged (ASP)",
              "and protonated (ASH), and a Na+ that binds the charged site but not the protonated",
              "one is the expected signature -- the pair of systems is only worth carrying",
              "forward if it produces a difference here or in the pocket.",
              ""]

    missing = [s for s, reps in systems.items() if len(reps) < 3]
    if missing:
        lines += [f"> **Note.** {', '.join(missing)} " +
                  ("has" if len(missing) == 1 else "have") +
                  " fewer than the 3 replicas the QC plan specifies, so the replicate spread is",
                  "> under-determined and R-hat is unavailable or unreliable for them.", ""]

    report = paths.PRODUCT / f"{PREFIX}_simulation_analysis_{args.date}.md"
    report.write_text("\n".join(lines) + "\n")

    summary_json = args.analysis_root / f"summary_{args.date}.json"
    summary_json.write_text(json.dumps({
        "date": args.date,
        "systems": {s: {"replicas": [r.replica for r in reps],
                        "total_ns": sum(float(r.summary.get("length_ns", 0.0)) for r in reps),
                        "checks": [vars(c) for c in all_checks[s]]}
                    for s, reps in systems.items()},
    }, indent=2) + "\n")

    print(f"Wrote:\n  {qc_csv}\n  {occ_csv}\n  {report}\n  {summary_json}")
    n_fail = sum(len(v) for v in failed.values())
    print(f"{len(results)} replicas, {len(systems)} systems, {total_ns:.0f} ns; "
          f"{n_fail} failed check(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

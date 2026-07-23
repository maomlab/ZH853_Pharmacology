#!/usr/bin/env python3
"""Predict receptor protonation states with PROPKA (Phase 2).

Runs PROPKA on MOR chain R, parses per-residue pKa, and recommends protonation states at
pH 7.4 -- flagging non-standard states, His tautomer sites, and the ambiguous D2.50 sodium
pocket. Writes product/02.02.00_protonation_*.md.

Run: ``python src/02.02.00_protonation.py``  (or ``make prep-protonation``).
"""

from __future__ import annotations

import subprocess
import sys
from datetime import date

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

from zh853mor import paths, prep, structure  # noqa: E402

PH = 7.4
IONIZABLE = {"ASP", "GLU", "HIS", "LYS", "ARG", "CYS", "TYR"}


def state_at_ph(resname: str, pka: float) -> tuple[str, bool]:
    """Return (recommended state, is_nonstandard) at pH 7.4."""
    if resname in ("ASP", "GLU"):
        return ("protonated/neutral (ASH/GLH)", True) if pka > PH else ("charged (-1)", False)
    if resname == "HIS":
        return ("HIP (+1)", True) if pka > PH else ("neutral (HID/HIE)", False)
    if resname == "LYS":
        return ("neutral (LYN)", True) if pka < PH else ("charged (+1)", False)
    if resname == "CYS":
        return ("thiolate (-1)", True) if pka < PH else ("thiol (neutral)", False)
    if resname == "ARG":
        return ("neutral", True) if pka < PH else ("charged (+1)", False)
    if resname == "TYR":
        return ("phenolate (-1)", True) if pka < PH else ("neutral", False)
    return ("standard", False)


def parse_pka(pka_file) -> list[tuple[str, int, float, float]]:
    """Parse the PROPKA SUMMARY section -> (resname, resid, pKa, model_pKa)."""
    rows, in_summary = [], False
    for line in pka_file.read_text().splitlines():
        if "SUMMARY OF THIS PREDICTION" in line:
            in_summary = True
            continue
        if in_summary:
            if line.strip().startswith("Free energy") or not line.strip():
                if rows:
                    break
                continue
            parts = line.split()
            if len(parts) >= 5 and parts[0] in IONIZABLE:
                try:
                    rows.append((parts[0], int(parts[1]), float(parts[3]), float(parts[4])))
                except ValueError:
                    continue
    return rows


def main() -> int:
    today = f"{date.today():%Y%m%d}"
    work = paths.ensure_dir(paths.INTERMEDIATE / "02.02.00_protonation")
    u = structure.load(paths.CRYOEM_PDB)
    rec_pdb = work / "receptorR.pdb"
    u.select_atoms("segid R and protein").write(str(rec_pdb))

    subprocess.run(["propka3", rec_pdb.name], cwd=str(work), check=True,
                   capture_output=True, text=True)
    rows = parse_pka(work / "receptorR.pka")

    disulf = {i for i, _, _ in structure.disulfides(u)} | {j for _, j, _ in structure.disulfides(u)}
    nonstd, his_sites = [], []
    for resname, resid, pka, _model in rows:
        state, ns = state_at_ph(resname, pka)
        if resname == "CYS" and resid in disulf:
            continue  # disulfide, not titratable
        if resname == "HIS":
            his_sites.append((resid, pka, state))
        if ns:
            nonstd.append((resname, resid, pka, state))

    lines = [f"# Receptor protonation (PROPKA, pH {PH}) — {today}", "",
             f"PROPKA on MOR chain R; {len(rows)} ionizable residues assessed. "
             "Default: standard states at pH 7.4 except those flagged below.", ""]

    lines.append("## Key functional / pocket residues")
    key = {prep.D250_SODIUM: "D2.50 sodium pocket", 149: "D3.32 anchor", 231: "E231 ECL2 (ZH853)",
           166: "D3.49 (DRY)", 299: "H6.52", 321: "H7.36 (ZH853-distinctive)"}
    by_id = {resid: (resname, pka) for resname, resid, pka, _ in rows}
    for resid, label in key.items():
        if resid in by_id:
            rn, pka = by_id[resid]
            state, _ = state_at_ph(rn, pka)
            lines.append(f"- **{rn}{resid}** ({label}): pKa {pka:.2f} → {state}")

    lines.append(f"\n## D2.50 (Asp116) — the load-bearing ambiguity")
    d250 = by_id.get(prep.D250_SODIUM)
    if d250:
        lines.append(f"pKa = {d250[1]:.2f}, i.e. ~at physiological pH. Occupancy of the proton is "
                     "genuinely uncertain → **build parallel systems (protonated ASH vs charged) and "
                     "compare** (SPECIFICATION D-11). Constant-pH MD is the rigorous fallback.")

    lines.append(f"\n## Non-standard protonation states to apply: {len(nonstd)}")
    if nonstd:
        lines.append("| Residue | pKa | recommended state |")
        lines.append("|---|---|---|")
        for rn, rid, pka, state in sorted(nonstd, key=lambda t: t[1]):
            lines.append(f"| {rn}{rid} | {pka:.2f} | {state} |")
    else:
        lines.append("None beyond D2.50 — all other ionizables take standard states.")

    lines.append(f"\n## Histidine tautomers (chain R): {len(his_sites)}")
    lines.append("PROPKA gives protonation, not tautomer; assign HID vs HIE from local H-bonding "
                 "(check the ligand/water network for H299 and H321).")
    for rid, pka, state in sorted(his_sites):
        lines.append(f"- HIS{rid}: pKa {pka:.2f} → {state}")

    (paths.PRODUCT / f"02.02.00_protonation_{today}.md").write_text("\n".join(lines) + "\n")
    print(f"Assessed {len(rows)} ionizables; {len(nonstd)} non-standard; D2.50 pKa="
          f"{by_id.get(prep.D250_SODIUM, ('', float('nan')))[1]:.2f}")
    print(f"Report: product/02.02.00_protonation_{today}.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

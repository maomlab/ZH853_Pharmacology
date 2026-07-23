#!/usr/bin/env python3
"""Comparative protein-ligand interaction fingerprints (Objective 1).

Computes heavy-atom interaction fingerprints for the ZH853 complex and every fetched
comparator, aligns them onto human OPRM1 numbering, and identifies which pocket contacts
are shared vs ZH853-distinctive. Writes a matrix (CSV), a narrative report (MD), a detailed
ZH853 fingerprint (MD), and a heatmap (PNG) to product/.

Run: ``python src/03.01.00_interaction_fingerprints.py``  (or ``make interactions``).
"""

from __future__ import annotations

import sys
from datetime import date

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from rdkit import RDLogger  # noqa: E402

from zh853mor import comparators, interactions, paths, structure  # noqa: E402

RDLogger.DisableLog("rdApp.*")

# Interaction -> (single-letter code, strength for the heatmap). Higher = stronger/more specific.
STRENGTH = {"ionic": 4.0, "cation_pi": 3.5, "aromatic": 3.0, "hbond": 2.0, "hydrophobic": 1.0}
CODE = {"ionic": "I", "cation_pi": "P", "aromatic": "A", "hbond": "H", "hydrophobic": "h"}

# Display order: ZH853 first, then peptide parents, other agonists, then antagonist.
ORDER = ["ZH853", "8F7R", "8EFQ", "6DDE", "8F7Q", "9WST", "9WSV",
         "5C1M", "8EF5", "8EFB", "8EFL", "8EFO", "7T2G", "4DKL"]


def cell_strength(fp: dict[int, interactions.ResidueFingerprint], hid: int) -> float:
    """Heatmap value for one residue in one structure."""
    r = fp.get(hid)
    if r is None:
        return 0.0
    if not r.interactions:
        return 0.5  # contact-only (within 4.5 A but no typed interaction)
    return max(STRENGTH[i] for i in r.interactions)


def cell_code(fp: dict[int, interactions.ResidueFingerprint], hid: int) -> str:
    r = fp.get(hid)
    if r is None:
        return ""
    if not r.interactions:
        return "."
    return "".join(CODE[i] for i in sorted(r.interactions, key=lambda x: -STRENGTH[x]))


def main() -> int:
    today = f"{date.today():%Y%m%d}"
    fps: dict[str, dict[int, interactions.ResidueFingerprint]] = {}
    names: dict[str, str] = {}
    classes: dict[str, str] = {}
    for pdb in ORDER:
        cx = comparators.load_complex(pdb)
        fps[pdb] = interactions.fingerprint(cx.receptor, cx.ligand, cx.offset_to_human)
        names[pdb] = cx.name
        classes[pdb] = cx.ligand_class
        print(f"  {pdb:6} {cx.name:18} {len(fps[pdb])} contact residues")

    # Union of pocket residues, sorted by TM/BW then residue id.
    resids = sorted({h for fp in fps.values() for h in fp})
    # canonical resname (from ZH853 structure preferred, else first structure that has it)
    resname: dict[int, str] = {}
    for hid in resids:
        for pdb in ORDER:
            if hid in fps[pdb]:
                resname[hid] = fps[pdb][hid].resname
                break

    agonists = [p for p in ORDER if p not in ("ZH853", "4DKL")]

    # ---- CSV matrix ----------------------------------------------------------
    paths.ensure_dir(paths.PRODUCT)
    csv = paths.PRODUCT / f"03.01.00_interaction_matrix_{today}.csv"
    header = ["resid", "resname", "BW"] + ORDER
    rows = [",".join(header)]
    for hid in resids:
        cells = [cell_code(fps[p], hid) for p in ORDER]
        rows.append(",".join([str(hid), resname[hid], structure.bw(hid), *cells]))
    csv.write_text("\n".join(rows) + "\n")

    # ---- shared vs ZH853-distinctive ----------------------------------------
    zh = fps["ZH853"]
    distinctive: list[tuple[int, int, str]] = []  # (n_agonists_sharing, hid, ztypes)
    for hid, r in zh.items():
        n_share = sum(hid in fps[p] for p in agonists)
        distinctive.append((n_share, hid, ",".join(sorted(r.interactions)) or "contact"))
    distinctive.sort(key=lambda t: (t[0], t[1]))

    lines = [f"# ZH853 vs MOR agonists — comparative interaction analysis ({today})", ""]
    lines.append(f"Heavy-atom geometric fingerprints (cutoffs: ionic {interactions.IONIC_CUT} A, "
                 f"H-bond {interactions.HBOND_CUT} A, hydrophobic {interactions.HYDROPHOBIC_CUT} A, "
                 f"pi-stack {interactions.ARO_CUT} A, cation-pi {interactions.CATIONPI_CUT} A). "
                 f"Codes: I=ionic P=cation-pi A=pi-stack H=hbond h=hydrophobic .=contact-only.")
    lines.append(f"\nAgonist reference set ({len(agonists)}): "
                 + ", ".join(f"{p}/{names[p]}" for p in agonists))
    lines.append(f"\n## ZH853 contacts ranked by distinctiveness (fewest agonists sharing first)\n")
    lines.append("| Residue | BW | ZH853 interactions | # agonists sharing | shared with |")
    lines.append("|---|---|---|---|---|")
    for n_share, hid, ztypes in distinctive:
        shared = [names[p] for p in agonists if hid in fps[p]]
        shared_s = ", ".join(shared) if shared else "**none (ZH853-unique)**"
        lines.append(f"| {resname[hid]}{hid} | {structure.bw(hid)} | {ztypes} | "
                     f"{n_share}/{len(agonists)} | {shared_s} |")

    report = paths.PRODUCT / f"03.01.00_interaction_comparison_{today}.md"
    report.write_text("\n".join(lines) + "\n")

    # ---- detailed ZH853 fingerprint -----------------------------------------
    zlines = [f"# ZH853 interaction fingerprint ({today})", "",
              "Receptor residues contacting ZH853 (human OPRM1 numbering).", ""]
    zlines.append("| Residue | BW | interactions | min dist (A) |")
    zlines.append("|---|---|---|---|")
    for hid in sorted(zh):
        r = zh[hid]
        types = ", ".join(sorted(r.interactions)) or "contact-only"
        zlines.append(f"| {r.resname}{hid} | {structure.bw(hid)} | {types} | {r.min_dist:.2f} |")
    zfile = paths.PRODUCT / f"03.01.00_zh853_fingerprint_{today}.md"
    zfile.write_text("\n".join(zlines) + "\n")

    # ---- heatmap -------------------------------------------------------------
    mat = np.array([[cell_strength(fps[p], hid) for p in ORDER] for hid in resids])
    fig, ax = plt.subplots(figsize=(0.55 * len(ORDER) + 3, 0.30 * len(resids) + 2))
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=4)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([f"{names[p]}\n({classes[p][:4]})" for p in ORDER],
                       rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(resids)))
    ax.set_yticklabels([f"{resname[h]}{h} {structure.bw(h)}" for h in resids], fontsize=6)
    for i, hid in enumerate(resids):
        for j, p in enumerate(ORDER):
            c = cell_code(fps[p], hid)
            if c and c != ".":
                ax.text(j, i, c, ha="center", va="center", fontsize=5,
                        color="black" if mat[i, j] < 3 else "white")
    ax.axvline(0.5, color="navy", lw=1.5)  # separate ZH853
    cbar = fig.colorbar(im, ax=ax, shrink=0.5)
    cbar.set_label("interaction strength", fontsize=7)
    ax.set_title("ZH853 vs MOR agonists: pocket interaction fingerprint", fontsize=9)
    fig.tight_layout()
    png = paths.PRODUCT / f"03.01.00_fingerprint_heatmap_{today}.png"
    fig.savefig(png, dpi=200)
    plt.close(fig)

    print(f"\nWrote:\n  {csv}\n  {report}\n  {zfile}\n  {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

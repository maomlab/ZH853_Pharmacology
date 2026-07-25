#!/usr/bin/env python3
"""Trim white borders from the MolStar renders for tight manuscript panels.

Crops each render to the bounding box of non-white content plus a small uniform margin,
and writes the trimmed image to product/manuscript/figures/.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
PRODUCT = HERE.parents[1] / "product"
FIGDIR = PRODUCT / "manuscript" / "figures"

# (source render, destination manuscript figure)
PAIRS = [
    ("03.10.00_molstar_overview_20260723.png", "fig6_molstar_overview.png"),
    ("03.10.00_molstar_pocket_20260723.png", "fig7_molstar_pocket.png"),
    ("03.10.00_molstar_membrane_20260725.png", "fig_membrane_placement.png"),
]


def trim(src: Path, dst: Path, threshold: int = 248, margin_frac: float = 0.02,
         inset_frac: float = 0.03) -> None:
    img = Image.open(src).convert("RGB")
    # drop the faint viewport border frame before content-trimming
    inset = int(inset_frac * min(img.width, img.height))
    img = img.crop((inset, inset, img.width - inset, img.height - inset))
    arr = np.asarray(img)
    # non-white = any channel below threshold
    mask = (arr < threshold).any(axis=2)
    ys, xs = np.where(mask)
    if len(xs) == 0:
        img.save(dst)
        return
    x0, x1, y0, y1 = xs.min(), xs.max(), ys.min(), ys.max()
    m = int(margin_frac * max(x1 - x0, y1 - y0))
    x0 = max(0, x0 - m)
    y0 = max(0, y0 - m)
    x1 = min(arr.shape[1] - 1, x1 + m)
    y1 = min(arr.shape[0] - 1, y1 + m)
    img.crop((x0, y0, x1 + 1, y1 + 1)).save(dst)
    print(f"{dst.name}: {img.width}x{img.height} -> {x1 - x0 + 1}x{y1 - y0 + 1} "
          f"(aspect {(x1 - x0 + 1) / (y1 - y0 + 1):.2f})")


def main() -> int:
    for s, d in PAIRS:
        trim(PRODUCT / s, FIGDIR / d)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

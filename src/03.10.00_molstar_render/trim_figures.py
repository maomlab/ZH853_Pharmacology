#!/usr/bin/env python3
"""Post-process the MolStar renders into tight manuscript figures.

- overview / pocket: crop to the non-white content bounding box (drop the viewport frame first).
- membrane: draw a solid bilayer band BEHIND the protein (spanning the cholesterol z-extent,
  full width) then crop -- this shows the inferred membrane explicitly as a sanity check and
  fills the lateral whitespace, while the molecule renders on top.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
PRODUCT = HERE.parents[1] / "product"
FIGDIR = PRODUCT / "manuscript" / "figures"

WHITE_THRESH = 248
INSET_FRAC = 0.03
MARGIN_FRAC = 0.02
BAND_RGB = (228, 210, 170)  # soft tan bilayer band


def _crop_to_content(arr: np.ndarray, margin_frac: float = MARGIN_FRAC) -> np.ndarray:
    mask = (arr < WHITE_THRESH).any(axis=2)
    ys, xs = np.where(mask)
    if len(xs) == 0:
        return arr
    x0, x1, y0, y1 = xs.min(), xs.max(), ys.min(), ys.max()
    m = int(margin_frac * max(x1 - x0, y1 - y0))
    x0, y0 = max(0, x0 - m), max(0, y0 - m)
    x1, y1 = min(arr.shape[1] - 1, x1 + m), min(arr.shape[0] - 1, y1 + m)
    return arr[y0:y1 + 1, x0:x1 + 1]


def _inset(im: Image.Image) -> Image.Image:
    d = int(INSET_FRAC * min(im.width, im.height))
    return im.crop((d, d, im.width - d, im.height - d))


def trim(src: Path, dst: Path) -> None:
    arr = np.asarray(_inset(Image.open(src).convert("RGB")))
    out = _crop_to_content(arr)
    Image.fromarray(out).save(dst)
    print(f"{dst.name}: trimmed -> {out.shape[1]}x{out.shape[0]}")


def composite_membrane(src: Path, dst: Path, margin_frac: float = 0.03) -> None:
    """Paint a bilayer band (cholesterol z-extent) behind the molecule, cropped TIGHT to the protein.

    The crop is the molecule's bounding box (+ a small margin), NOT the band's -- so the band only
    fills the tight crop and does not leave wide membrane-only strips on the sides.
    """
    arr = np.asarray(_inset(Image.open(src).convert("RGB"))).copy()
    molecule = (arr < 246).any(axis=2)                                     # protein + cholesterol + ligand
    ys, xs = np.where(molecule)
    if len(xs) == 0:
        Image.fromarray(arr).save(dst)
        return
    x0, x1, y0, y1 = xs.min(), xs.max(), ys.min(), ys.max()               # crop box = the molecule
    mx, my = int(margin_frac * (x1 - x0)), int(margin_frac * (y1 - y0))
    x0, y0 = max(0, x0 - mx), max(0, y0 - my)
    x1, y1 = min(arr.shape[1] - 1, x1 + mx), min(arr.shape[0] - 1, y1 + my)

    r, g, b = arr[..., 0].astype(int), arr[..., 1].astype(int), arr[..., 2].astype(int)
    orange = (r > 165) & (g > 65) & (g < 195) & (b < 100) & (r - b > 85)   # cholesterol -> band rows
    rows = np.where(orange.any(axis=1))[0]
    if len(rows):
        band = np.zeros(arr.shape[:2], bool)
        band[int(rows.min()):int(rows.max()) + 1, :] = True
        arr[band & ~molecule] = BAND_RGB                                  # band behind the molecule

    out = arr[y0:y1 + 1, x0:x1 + 1]                                       # tight crop to the protein
    Image.fromarray(out).save(dst)
    print(f"{dst.name}: membrane band + tight crop -> {out.shape[1]}x{out.shape[0]}")


def main() -> int:
    trim(PRODUCT / "03.10.00_molstar_overview_20260723.png", FIGDIR / "fig6_molstar_overview.png")
    trim(PRODUCT / "03.10.00_molstar_pocket_20260723.png", FIGDIR / "fig7_molstar_pocket.png")
    composite_membrane(PRODUCT / "03.10.00_molstar_membrane_20260725.png",
                       FIGDIR / "fig_membrane_placement.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

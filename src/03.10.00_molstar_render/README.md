# MolStar headless 3D rendering

Renders the MOR–ZH853 complex to publication PNGs using the prebuilt **Mol\*** viewer driven by
headless Chrome (software WebGL) — no display or GPU required.

## Files
- `package.json` — pins `molstar` + `puppeteer` (installs a Chromium build).
- `viewer.html` — loads the Mol\* viewer bundle; white background, outline + ambient occlusion,
  UI/axis hidden; exposes `loadPDB`, `loadMVS`, and scene helpers.
- `build_pocket_mvs.py` — builds `pocket.mvsj` (a MolViewSpec scene): receptor cartoon (faded),
  ZH853 + key pocket residues as sticks with BW labels, and salt-bridge/H-bond distance lines.
- `render.js` — Puppeteer script: renders Scene 1 (whole complex overview, default preset) and
  Scene 2 (pocket + interactions from the MVS), screenshots each `#app` canvas.

## Run
```bash
npm install                     # one-time (downloads Chromium)
python build_pocket_mvs.py      # -> pocket.mvsj
node render.js                  # -> product/03.10.00_molstar_{overview,pocket}_*.png
```
or `make molstar-render` from the repo root. Outputs are copied to
`product/manuscript/figures/fig6_molstar_overview.png` and `fig7_molstar_pocket.png`.

## Notes
- `node_modules/` and `pocket.mvsj` are git-ignored (regenerate via the commands above).
- Chrome is launched with `--use-gl=angle --use-angle=swiftshader` for headless WebGL2; kill any
  orphaned `Chrome for Testing` processes if a run hangs, then re-run.
- Requires Node 18+.

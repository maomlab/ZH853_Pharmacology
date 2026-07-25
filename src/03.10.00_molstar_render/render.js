// Headless MolStar rendering of the MOR-ZH853 complex for the manuscript.
// Runs the prebuilt MolStar viewer in headless Chrome (software WebGL) and screenshots scenes.
//   node render.js
const fs = require('fs');
const path = require('path');
const puppeteer = require('puppeteer');

const HERE = __dirname;
const REPO = path.resolve(HERE, '../..');
const FULL_PDB = path.join(REPO, 'data', 'mor_gi_scfv16_refine-coot-22_real_space_refined_169_edit.pdb');
const RECLIG_PDB = path.join(REPO, 'intermediate', '02.01.00_components', 'receptor_ligand.pdb');
const OUTDIR = path.join(REPO, 'product');
const FIGDIR = path.join(OUTDIR, 'manuscript', 'figures');

async function shoot(page, name) {
  const el = await page.$('#app');
  const out = path.join(OUTDIR, name);
  await el.screenshot({ path: out });
  console.log('wrote', out);
  return out;
}

(async () => {
  const browser = await puppeteer.launch({
    headless: 'new',
    args: ['--no-sandbox', '--enable-webgl', '--ignore-gpu-blocklist',
      '--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader',
      '--in-process-gpu', '--window-size=1650,1250'],
  });
  const page = await browser.newPage();
  await page.setViewport({ width: 1650, height: 1250, deviceScaleFactor: 2 });
  page.on('console', m => { if (m.type() === 'error') console.log('PAGE ERR:', m.text()); });
  page.on('pageerror', e => console.log('PAGEERROR:', e.message));

  await page.goto('file://' + path.join(HERE, 'viewer.html'), { waitUntil: 'networkidle0' });
  await page.evaluate(() => initViewer());
  await page.waitForFunction('window.ready === true', { timeout: 60000 });
  console.log('viewer ready');

  // --- Scene 1: overview (full complex) ---
  const full = fs.readFileSync(FULL_PDB, 'utf8');
  await page.evaluate((s) => loadPDB(s), full);
  await page.waitForFunction('window.loaded === true', { timeout: 120000 });
  await page.evaluate(() => sceneOverview());
  const cam = JSON.parse(fs.readFileSync(path.join(HERE, 'overview_camera.json'), 'utf8'));
  await page.evaluate((c) => setOverviewCamera(c), cam);
  await shoot(page, '03.10.00_molstar_overview_20260723.png');

  // --- Scene 2: orthosteric pocket + interactions (MolViewSpec) ---
  await page.evaluate(() => { window.loaded = false; window.viewer.plugin.clear(); });
  const mvs = fs.readFileSync(path.join(HERE, 'pocket.mvsj'), 'utf8');
  await page.evaluate((s) => loadMVS(s), mvs);
  await page.waitForFunction('window.loaded === true', { timeout: 120000 });
  await page.evaluate(() => settle());
  const pocket = await shoot(page, '03.10.00_molstar_pocket_20260723.png');

  // --- Scene 3: membrane placement (MolViewSpec) ---
  let membrane = null;
  const memPath = path.join(HERE, 'membrane.mvsj');
  if (fs.existsSync(memPath)) {
    await page.evaluate(() => { window.loaded = false; window.viewer.plugin.clear(); });
    await page.evaluate((s) => loadMVS(s), fs.readFileSync(memPath, 'utf8'));
    await page.waitForFunction('window.loaded === true', { timeout: 120000 });
    await page.evaluate(() => settle());
    membrane = await shoot(page, '03.10.00_molstar_membrane_20260725.png');
  }

  await browser.close();
  // NOTE: manuscript figures are written by trim_figures.py (trim/composite), not here,
  // so the staged figures are always the processed versions.
  console.log('done');
})().catch(e => { console.error(e); process.exit(1); });

// Cartoon-quality movies of the production trajectories, rendered headlessly with MolStar.
//
// **Runs locally**, on the movie files 01_export_movie_trajectory.py wrote on the cluster. The
// companion to 02_render_dashboard.py: that one answers "what is this trajectory doing, and does
// the number agree?", this one answers "what does it LOOK like" -- cartoon helices, the ligand in
// sticks, the annular lipids around it, at a quality that can go in a talk.
//
//   npm install                 # once
//   node 03_render_molstar.js                      # every exported replica
//   node 03_render_molstar.js --system apo_ASH --replica prod_r1
//   node 03_render_molstar.js --fps 25 --width 1600 --height 900
//
// Frames are screenshotted one at a time and piped to ffmpeg. MolStar's own mp4 export is not
// used: it drives its own animation clock, so the frames it emits are not the frames the movie
// file actually holds, and the caption could not be kept in step with them.

const fs = require('fs');
const path = require('path');
const os = require('os');
const { spawnSync } = require('child_process');
const puppeteer = require('puppeteer');

const HERE = __dirname;
const REPO = path.resolve(HERE, '../..');
const MOVIE_ROOT = path.join(REPO, 'intermediate', '02.12.00_render_trajectory_movies');
const PRODUCT = path.join(REPO, 'product');
const PREFIX = '02.12.00';

const USAGE = `Render exported movie trajectories as MolStar cartoon movies.

  node 03_render_molstar.js [options]

  --system <name>     only this system, e.g. apo_ASH
  --replica <name>    only this replica, e.g. prod_r1
  --fps <n>           frame rate of the mp4            (default 25)
  --width <px>        canvas width, must be even       (default 1280)
  --height <px>       canvas height, must be even      (default 720)
  --zoom <f>          >1 pulls the camera back         (default 1.15)
  --movies <dir>      exported movie trajectories      (default intermediate/02.12.00_*)
  --out <dir>         where the mp4s go                (default product/)
`;

function parseArgs(argv) {
  const out = { fps: 25, width: 1280, height: 720, zoom: 1.15, movies: MOVIE_ROOT, out: PRODUCT };
  if (argv.includes('--help') || argv.includes('-h')) {
    process.stdout.write(USAGE);
    process.exit(0);
  }
  for (let i = 2; i < argv.length; i += 2) {
    const key = argv[i].replace(/^--/, '');
    const val = argv[i + 1];
    if (val === undefined) throw new Error(`--${key} needs a value`);
    out[key] = /^(fps|width|height|zoom)$/.test(key) ? Number(val) : val;
  }
  if (out.width % 2 || out.height % 2) {
    throw new Error('--width and --height must both be even (yuv420p cannot encode odd sizes)');
  }
  return out;
}

function findMovies(root, system, replica) {
  if (!fs.existsSync(root)) return [];
  const found = [];
  for (const dir of fs.readdirSync(root).sort()) {
    const sub = path.join(root, dir);
    if (!fs.statSync(sub).isDirectory()) continue;
    if (system && dir !== system) continue;
    for (const f of fs.readdirSync(sub).sort()) {
      if (!f.endsWith('_movie.json')) continue;
      if (replica && !f.startsWith(`${replica}_movie`)) continue;
      found.push(path.join(sub, f));
    }
  }
  return found;
}

function stamp() {
  const d = new Date();
  const p = n => String(n).padStart(2, '0');
  return `${d.getFullYear()}${p(d.getMonth() + 1)}${p(d.getDate())}`;
}

function haveFfmpeg() {
  return spawnSync('ffmpeg', ['-version'], { stdio: 'ignore' }).status === 0;
}

// One palette, shared with 02_render_dashboard.py so the two movies of the same replica are the
// same colours. Keys match what buildRepresentations() expects.
const COLORS = {
  receptor: 0x2b5d8a, other: 0x9aa3ab, ligand: 0xe8710a,
  lipid: 0xc9bda6, water: 0x3fa7c4, ion: 0x7d3c98,
};
const LEGEND_TEXT = {
  receptor: 'receptor (cartoon)', other: 'Gi / other chains', ligand: 'ligand',
  lipid: 'annular lipids + phosphates', water: 'pocket water', ion: 'ions',
};

function legendFor(drawn) {
  // Only what was actually drawn: an apo system has no ligand row, a receptor-only build no
  // lipid row, and a legend naming a colour that is not on screen is worse than none.
  return drawn.map(k => `<span style="background:#${COLORS[k].toString(16).padStart(6, '0')}">` +
    `</span>${LEGEND_TEXT[k]}`).join('<br>');
}

// Bounding sphere of the receptor and ligand chains, straight from the movie PDB's ATOM records.
function focusSphere(pdbText, chains) {
  const wanted = new Set([...(chains.receptor || []), chains.ligand].filter(Boolean));
  const pts = [];
  for (const line of pdbText.split('\n')) {
    if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) continue;
    if (!wanted.has(line.substr(21, 1))) continue;
    pts.push([+line.substr(30, 8), +line.substr(38, 8), +line.substr(46, 8)]);
  }
  if (!pts.length) throw new Error('no receptor or ligand atoms found in the movie PDB');
  const c = [0, 1, 2].map(i => pts.reduce((s, p) => s + p[i], 0) / pts.length);
  const radius = Math.max(...pts.map(
    p => Math.hypot(p[0] - c[0], p[1] - c[1], p[2] - c[2])));
  return { center: c, radius };
}

async function renderOne(page, jsonPath, opts) {
  const meta = JSON.parse(fs.readFileSync(jsonPath, 'utf8'));
  const dir = path.dirname(jsonPath);
  const pdb = fs.readFileSync(path.join(dir, meta.files.topology), 'utf8');
  const xtc = fs.readFileSync(path.join(dir, meta.files.trajectory)).toString('base64');
  const times = meta.movie.time_ns;

  await page.evaluate(() => { window.loaded = false; window.viewer.plugin.clear(); });
  const frameCount = await page.evaluate((p, x) => loadMovie(p, x), pdb, xtc);
  await page.waitForFunction('window.loaded === true', { timeout: 180000 });
  if (!frameCount) throw new Error('MolStar loaded no trajectory frames');
  if (frameCount !== times.length) {
    // Not fatal -- the caption just has to follow the file that was actually loaded.
    console.log(`  note: XTC holds ${frameCount} frames, the JSON lists ${times.length}`);
  }
  const chains = (meta.atoms && meta.atoms.chains) || {};
  const shown = await page.evaluate((c, col) => buildRepresentations(c, col), chains, COLORS);
  // Camera ONCE, before the first frame, and never touched again: a per-frame reset would
  // follow the receptor's own drift and hide exactly the motion the movie is being made for.
  await page.evaluate((f, z) => setMembraneCamera(f, z), focusSphere(pdb, chains), opts.zoom);
  await page.evaluate(l => setLegend(l), legendFor(shown));

  const outName = `${PREFIX}_${meta.build}_${meta.replica}_molstar_${stamp()}.mp4`;
  const outPath = path.join(opts.out, outName);
  const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'zh853-molstar-'));
  const el = await page.$('#shot');
  const aligned = meta.movie.aligned
    ? `superposed on ${meta.movie.align_selection}`
    : 'NOT superposed — box drift and tilt left in';

  const n = Math.min(frameCount, times.length);
  for (let i = 0; i < n; i++) {
    await page.evaluate(k => setFrame(k), i);
    await page.evaluate(html => setCaption(html),
      `<b>${meta.build} / ${meta.replica}</b>   t = ${times[i].toFixed(1)} ns   ` +
      `frame ${i + 1}/${n}\n${aligned}`);
    await page.evaluate(() => drawn());
    await el.screenshot({ path: path.join(tmp, `f${String(i).padStart(5, '0')}.png`) });
    if ((i + 1) % 50 === 0 || i + 1 === n) process.stdout.write(`\r  ${i + 1}/${n} frames`);
  }
  process.stdout.write('\n');

  const ff = spawnSync('ffmpeg', [
    '-y', '-loglevel', 'error', '-framerate', String(opts.fps),
    '-i', path.join(tmp, 'f%05d.png'),
    '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2',
    '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '20', '-preset', 'medium', outPath,
  ], { stdio: 'inherit' });
  fs.rmSync(tmp, { recursive: true, force: true });
  if (ff.status !== 0) throw new Error(`ffmpeg exited ${ff.status}`);
  console.log(`  wrote ${outPath}  [${n} frames, ${(n / opts.fps).toFixed(1)} s]`);
  return outPath;
}

(async () => {
  const opts = parseArgs(process.argv);
  if (!haveFfmpeg()) {
    console.error('ERROR: ffmpeg is not on PATH. `conda install -c conda-forge ffmpeg`.');
    process.exit(1);
  }
  const movies = findMovies(opts.movies, opts.system, opts.replica);
  if (!movies.length) {
    console.error(`ERROR: no exported movies under ${opts.movies}. Run ` +
      '01_export_movie_trajectory.py on the cluster and copy the directory back.');
    process.exit(1);
  }
  fs.mkdirSync(opts.out, { recursive: true });

  const browser = await puppeteer.launch({
    headless: 'new',
    args: ['--no-sandbox', '--enable-webgl', '--ignore-gpu-blocklist',
      '--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader',
      '--in-process-gpu', `--window-size=${opts.width},${opts.height}`],
  });
  const page = await browser.newPage();
  await page.setViewport({ width: opts.width, height: opts.height, deviceScaleFactor: 1 });
  page.on('pageerror', e => console.log('PAGEERROR:', e.message));
  page.on('console', m => { if (m.type() === 'error') console.log('PAGE ERR:', m.text()); });

  await page.goto('file://' + path.join(HERE, 'movie_viewer.html'), { waitUntil: 'networkidle0' });
  await page.evaluate((w, h) => initViewer(w, h), opts.width, opts.height);
  await page.waitForFunction('window.ready === true', { timeout: 120000 });

  let failures = 0;
  for (const js of movies) {
    console.log(`--- ${path.basename(path.dirname(js))} / ${path.basename(js)}`);
    try {
      await renderOne(page, js, opts);
    } catch (err) {
      // One unusable replica must not take the rest of the panel with it.
      console.error(`FAILED ${js}: ${err.message}`);
      failures += 1;
    }
  }
  await browser.close();
  process.exit(failures ? 1 : 0);
})();

#!/usr/bin/env python3
"""Render a movie trajectory as a diagnostic dashboard: structure beside its own time series.

**Runs locally** (`zh853mor-local`) on the few MB per replica that `01_export_movie_trajectory.py`
copied back. No trajectory, no cluster.

The point of this figure is not to look like the receptor. It is to put the 3D motion and the
scalar traces on the SAME clock, because that is the join a static figure cannot make: a step in
the ligand RMSD at 230 ns is a number until you can see, in the same instant, that the macrocycle
has rotated out past W6.48 and two annular lipids have moved into the vacated shell.

    +---------------+------------+-----------+
    |               |  pocket    |  TM RMSD  |
    |  membrane     |  (x-y)     +-----------+
    |  (x-z)        |            |  lig RMSD |
    |  receptor in  +------------+-----------+
    |  the bilayer  |  contacts  |  bilayer  |
    |               |            |  Na+/2.50 |
    +---------------+------------+-----------+

Every panel is derived from the movie files themselves, so the dashboard renders for a replica
that has not been through `02.11.00` at all. Where the reduced `.npz` from that stage IS present
it is preferred for the traces, because those were measured on every frame of the full-resolution
trajectory rather than on the ~300 the movie kept, and the ligand RMSD there is against the
DEPOSITED pose rather than against frame 0. The panel titles say which was used.

    python 02_render_dashboard.py                      # every exported replica
    python 02_render_dashboard.py --system apo_ASH --replica prod_r1
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import MDAnalysis as mda  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import animation  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

from zh853mor import md, paths  # noqa: E402

MOVIE_ROOT = paths.INTERMEDIATE / "02.12.00_render_trajectory_movies"
ANALYSIS_ROOT = paths.INTERMEDIATE / "02.11.00_analyze_simulations"
PREFIX = "02.12.00"

FPS = 25
DPI = 100
FIGSIZE = (16.0, 9.0)          # 1600 x 900 at DPI 100
CONTACT_CUT = 4.5              # A, heavy-atom: the same criterion as 03.01.00 and 02.11.00
POCKET_PAD = 14.0              # A around the ligand for the pocket view
POCKET_SLAB = 9.0              # A in z: without it the whole bundle overlaps the pocket

C = {
    "receptor": "#2b5d8a", "other": "#9aa3ab", "ligand": "#e8710a",
    "lipid": "#cfc7b8", "phosphate": "#8a8f98", "water": "#3fa7c4", "ion": "#7d3c98",
    "contact": "#cc2f26", "idle": "#9aa3ab", "cursor": "#cc2f26", "grid": "#dddddd",
}


# --- loading -------------------------------------------------------------------------------------
class Movie:
    """One exported replica: coordinates in memory, plus everything the panels need.

    The whole movie is held as one (frames, atoms, 3) array. At the default 300 frames and ~13k
    atoms that is ~47 MB -- cheaper than re-reading the XTC 300 times, and it lets every panel's
    axis limits be set from the FULL extent before the first frame is drawn, so nothing jumps.
    """

    def __init__(self, json_path: Path) -> None:
        self.meta = json.loads(json_path.read_text())
        d = json_path.parent
        self.u = mda.Universe(str(d / self.meta["files"]["topology"]),
                              str(d / self.meta["files"]["trajectory"]))
        self.system = str(self.meta["build"])
        self.replica = str(self.meta["replica"])
        self.time_ns = np.asarray(self.meta["movie"]["time_ns"], dtype=float)
        g = self.meta["atoms"]["groups"]
        self.slices = {k: slice(v[0], v[1]) for k, v in g.items()}
        self.xyz = np.stack([self.u.atoms.positions.copy() for _ in self.u.trajectory])
        self.protein = self.u.atoms[self.slices["protein"]]
        self.ligand = self.u.atoms[self.slices["ligand"]]

        # Chains come from the exported PDB, where each polypeptide was written as its own. Which
        # of them is the RECEPTOR was settled by the exporter and travels in the JSON, so this
        # figure and the MolStar render cannot disagree about what the subject of the movie is.
        self.ca = self.protein.select_atoms("name CA")
        self.ca_ix = self.ca.ix
        self.chain_of_ca = np.array([a.chainID for a in self.ca])
        named = self.meta["atoms"].get("chains", {}).get("receptor")
        self.receptor_chains = set(named) if named else set(self.chain_of_ca[:1])

        # Key residues (pocket shell + the functional positions), one heavy-atom group each.
        self.key: dict[int, np.ndarray] = {}
        for resid in md.KEY_RESIDUES:
            sel = self.protein.select_atoms(f"resid {resid}")
            if len(sel):
                self.key[resid] = sel.ix
        self.key_labels = {r: md.KEY_RESIDUES[r] for r in self.key}

    @property
    def n_frames(self) -> int:
        return self.xyz.shape[0]

    def contacts(self) -> tuple[np.ndarray, list[int]]:
        """(frames x residues) boolean contact with the ligand, at the project's 4.5 A criterion."""
        resids = sorted(self.key)
        if not len(self.ligand) or not resids:
            return np.zeros((self.n_frames, 0), dtype=bool), []
        lig = self.xyz[:, self.ligand.ix, :]
        out = np.zeros((self.n_frames, len(resids)), dtype=bool)
        for j, r in enumerate(resids):
            res = self.xyz[:, self.key[r], :]
            # (frames, res_atoms, lig_atoms) -- small enough to do in one shot per residue
            d = np.linalg.norm(res[:, :, None, :] - lig[:, None, :, :], axis=-1)
            out[:, j] = d.min(axis=(1, 2)) < CONTACT_CUT
        return out, resids

    def thickness(self) -> np.ndarray:
        """Phosphate-to-phosphate bilayer thickness per frame, from the movie's own headgroups."""
        ix = np.r_[self.u.atoms[self.slices["phosphate"]].ix,
                   self.u.atoms[self.slices["lipid"]].select_atoms("name P").ix]
        if ix.size < 10:
            return np.full(self.n_frames, np.nan)
        return np.array([md.bilayer_thickness(self.xyz[i, ix, :]) for i in range(self.n_frames)])

    def rmsd_ca(self) -> np.ndarray:
        ref = self.xyz[0, self.ca_ix, :]
        d = self.xyz[:, self.ca_ix, :] - ref
        return np.sqrt((d ** 2).sum(axis=2).mean(axis=1))

    def rmsd_ligand(self) -> np.ndarray:
        if not len(self.ligand):
            return np.full(self.n_frames, np.nan)
        ref = self.xyz[0, self.ligand.ix, :]
        d = self.xyz[:, self.ligand.ix, :] - ref
        return np.sqrt((d ** 2).sum(axis=2).mean(axis=1))


def reduced_series(system: str, replica: str, root: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Full-resolution traces from `02.11.00` for this replica, or {} if it has not been reduced."""
    js = root / system / f"{replica}.json"
    npz = js.with_suffix(".npz")
    if not (js.exists() and npz.exists()):
        return {}
    try:
        result = md.ReplicaResult(json.loads(js.read_text()), npz)
    except ValueError:
        return {}
    out = {}
    for name in ("rmsd_ca_tm", "lig_rmsd_pose", "thickness", "na_d250_dist"):
        try:
            t, y = result.timeseries(name)
        except KeyError:
            continue
        if np.isfinite(y).any():
            out[name] = (t, y)
    out["t0_ns"] = (np.array([result.t0_ns]), np.array([np.nan]))  # marker, not a series
    return out


# --- drawing -------------------------------------------------------------------------------------
def ca_segments(pos: np.ndarray, chains: np.ndarray, ax0: int, ax1: int) -> np.ndarray:
    """Consecutive CA-CA segments projected onto two axes, never crossing a chain break.

    Drawing the trace as one polyline would join the receptor's C-terminus to Gai's N-terminus
    with a line straight across the picture, which is exactly the artefact the exported PDB's
    chain IDs exist to prevent.
    """
    same = chains[:-1] == chains[1:]
    a = pos[:-1][same][:, [ax0, ax1]]
    b = pos[1:][same][:, [ax0, ax1]]
    return np.stack([a, b], axis=1)


def style(ax, xlabel: str, ylabel: str, title: str = "") -> None:
    ax.set_xlabel(xlabel, fontsize=7)
    ax.set_ylabel(ylabel, fontsize=7)
    if title:
        ax.set_title(title, fontsize=8, loc="left")
    ax.tick_params(labelsize=6, length=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def series_panel(ax, movie: Movie, reduced, key: str, fallback: np.ndarray,
                 ylabel: str, reduced_title: str, own_title: str):
    """Plot one trace, preferring the reduced full-resolution version, and return its cursor."""
    if key in reduced:
        t, y = reduced[key]
        title = reduced_title
    else:
        t, y = movie.time_ns, fallback
        title = own_title
    finite = np.isfinite(y)
    if not finite.any():
        # `title` here is the fallback's own title, which says WHY there is nothing to draw --
        # an empty panel labelled only with its units leaves the reader guessing.
        ax.text(0.5, 0.5, title, ha="center", va="center", fontsize=7, style="italic",
                wrap=True, transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        return None
    ax.plot(t, y, lw=0.7, color=C["receptor"])
    ax.set_xlim(float(min(t.min(), movie.time_ns.min())),
                float(max(t.max(), movie.time_ns.max())))
    lo, hi = np.nanmin(y[finite]), np.nanmax(y[finite])
    pad = 0.08 * (hi - lo) + 1e-6
    ax.set_ylim(lo - pad, hi + pad)
    if "t0_ns" in reduced and reduced["t0_ns"][0][0] > 0:
        # Where 02.11.00 judged the relaxation to end; means are taken after it.
        ax.axvline(float(reduced["t0_ns"][0][0]), color="#888888", lw=0.6, ls=":")
    style(ax, "time (ns)", ylabel, title)
    return ax.axvline(movie.time_ns[0], color=C["cursor"], lw=1.1)


def build_figure(movie: Movie, reduced: dict):
    """Lay the dashboard out once and return (figure, per-frame update function)."""
    contacts, contact_resids = movie.contacts()
    thickness = movie.thickness()
    rmsd_ca, rmsd_lig = movie.rmsd_ca(), movie.rmsd_ligand()
    has_ligand = bool(len(movie.ligand))

    fig = plt.figure(figsize=FIGSIZE, dpi=DPI)
    gs = fig.add_gridspec(3, 3, width_ratios=[1.15, 1.15, 1.0], height_ratios=[1, 1, 1],
                          left=0.035, right=0.985, top=0.905, bottom=0.065,
                          wspace=0.24, hspace=0.42)
    ax_mem = fig.add_subplot(gs[:, 0])
    ax_pocket = fig.add_subplot(gs[0:2, 1])
    ax_strip = fig.add_subplot(gs[2, 1])
    ax_r1 = fig.add_subplot(gs[0, 2])
    ax_r2 = fig.add_subplot(gs[1, 2])
    # The bottom-right cell carries two short traces rather than one: bilayer thickness is the
    # membrane QC signal and is computable from the movie alone, while the Na+/D2.50 distance is
    # the direct test of the ASP/ASH pair the panel was built to compare (D-11) and only exists
    # once the replica has been reduced. Neither should displace the other.
    gs_bottom = gs[2, 2].subgridspec(2, 1, hspace=0.85)
    ax_r3 = fig.add_subplot(gs_bottom[0])
    ax_r4 = fig.add_subplot(gs_bottom[1])

    # --- membrane view: x-z, the whole system against the bilayer ------------------------------
    is_rec = np.isin(movie.chain_of_ca, list(movie.receptor_chains))
    lipid_ix = movie.u.atoms[movie.slices["lipid"]].ix
    phos_ix = movie.u.atoms[movie.slices["phosphate"]].ix
    water_ix = movie.u.atoms[movie.slices["water"]].select_atoms("name O").ix
    ion_ix = movie.u.atoms[movie.slices["ion"]].ix
    lig_ix = movie.ligand.ix

    lo = movie.xyz.reshape(-1, 3).min(axis=0)
    hi = movie.xyz.reshape(-1, 3).max(axis=0)
    ax_mem.set_xlim(lo[0] - 2, hi[0] + 2)
    ax_mem.set_ylim(lo[2] - 2, hi[2] + 2)
    ax_mem.set_aspect("equal")
    style(ax_mem, "x (A)", "z (A)  -- membrane normal",
          "membrane view: CA coloured by displacement from frame 0")

    mem_lipid = ax_mem.scatter([], [], s=1.0, c=C["lipid"], linewidths=0, zorder=1)
    mem_phos = ax_mem.scatter([], [], s=5.0, c=C["phosphate"], linewidths=0, zorder=2)
    mem_water = ax_mem.scatter([], [], s=6.0, c=C["water"], linewidths=0, zorder=3)
    mem_ion = ax_mem.scatter([], [], s=14.0, c=C["ion"], linewidths=0, zorder=4)
    # The receptor trace carries the diagnostic colour: per-residue displacement from frame 0
    # says WHICH part of the bundle is moving, which no scalar RMSD can. Grey -> red rather than
    # a stock colormap, because both ends have to read against a white page: a pale yellow low
    # end would make the parts that are behaving invisible instead of quiet.
    disp_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=5.0)
    disp_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "displacement", ["#8d99a6", "#f0ad4e", "#cc2f26", "#54120f"])
    mem_rec = LineCollection([], linewidths=1.6, cmap=disp_cmap, norm=disp_norm, zorder=6)
    mem_oth = LineCollection([], linewidths=1.0, colors=C["other"], zorder=5)
    ax_mem.add_collection(mem_rec)
    ax_mem.add_collection(mem_oth)
    mem_lig = LineCollection([], linewidths=2.2, colors=C["ligand"], zorder=7)
    ax_mem.add_collection(mem_lig)
    mem_lig_pts = ax_mem.scatter([], [], s=9.0, c=C["ligand"], linewidths=0, zorder=8)
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=disp_norm, cmap=disp_cmap),
                      ax=ax_mem, fraction=0.03, pad=0.01)
    cb.set_label("CA displacement (A)", fontsize=6)
    cb.ax.tick_params(labelsize=5)

    # --- pocket view: x-y, zoomed on the ligand -------------------------------------------------
    centre_ix = lig_ix if has_ligand else np.concatenate(
        [movie.key[r] for r in sorted(movie.key)]) if movie.key else movie.ca_ix
    centre = movie.xyz[:, centre_ix, :].mean(axis=1).mean(axis=0)
    ax_pocket.set_xlim(centre[0] - POCKET_PAD, centre[0] + POCKET_PAD)
    ax_pocket.set_ylim(centre[1] - POCKET_PAD, centre[1] + POCKET_PAD)
    ax_pocket.set_aspect("equal")
    style(ax_pocket, "x (A)", "y (A)",
          "pocket, extracellular view -- red = within 4.5 A of the ligand")
    pk_lipid = ax_pocket.scatter([], [], s=4.0, c=C["lipid"], linewidths=0, zorder=1)
    pk_water = ax_pocket.scatter([], [], s=18.0, c=C["water"], linewidths=0, zorder=2)
    pk_ion = ax_pocket.scatter([], [], s=34.0, c=C["ion"], linewidths=0, zorder=3)
    pk_trace = LineCollection([], linewidths=1.0, colors=C["other"], zorder=4)
    ax_pocket.add_collection(pk_trace)
    # facecolors=, not c=: `c` installs an (empty) scalar array which the draw pass re-maps
    # through the colormap, silently discarding every per-frame set_facecolor and leaving the
    # residue markers invisible.
    pk_key = ax_pocket.scatter([], [], s=42.0, facecolors=C["idle"], linewidths=0.4,
                               edgecolors="white", zorder=6)
    pk_lig = LineCollection([], linewidths=2.4, colors=C["ligand"], zorder=7)
    ax_pocket.add_collection(pk_lig)
    pk_labels = [ax_pocket.text(0, 0, "", fontsize=5.5, ha="center", va="bottom",
                                zorder=8, color="#333333") for _ in movie.key]
    lig_bonds = np.array([[b[0].ix, b[1].ix] for b in movie.ligand.bonds]) \
        if has_ligand and len(movie.ligand.bonds) else np.empty((0, 2), dtype=int)

    # --- contact strip: the occupancy analysis, accumulating as the movie plays ------------------
    if contact_resids:
        # Two flat colours, not a ramp: contact is a yes/no at 4.5 A, and a sequential colormap
        # tints the empty background pink, which reads as "weakly in contact everywhere".
        binary = matplotlib.colors.ListedColormap(["#ffffff", C["contact"]])
        ax_strip.imshow(contacts.T, aspect="auto", origin="lower", cmap=binary,
                        vmin=0, vmax=1, interpolation="nearest",
                        extent=(movie.time_ns[0], movie.time_ns[-1], -0.5,
                                len(contact_resids) - 0.5))
        ax_strip.set_yticks(range(len(contact_resids)))
        ax_strip.set_yticklabels([f"{r} {movie.key_labels[r]}" for r in contact_resids],
                                 fontsize=4.5)
        ax_strip.set_xlim(float(movie.time_ns[0]), float(movie.time_ns[-1]))
        for y in np.arange(len(contact_resids)) - 0.5:
            ax_strip.axhline(y, color=C["grid"], lw=0.3)
        strip_cursor = ax_strip.axvline(movie.time_ns[0], color="#222222", lw=1.1)
        style(ax_strip, "time (ns)", "", "ligand contacts (movie frames, 4.5 A)")
    else:
        ax_strip.text(0.5, 0.5, "apo system: no ligand contacts", ha="center", va="center",
                      fontsize=7, style="italic", transform=ax_strip.transAxes)
        ax_strip.set_xticks([])
        ax_strip.set_yticks([])
        strip_cursor = None

    # --- the traces ------------------------------------------------------------------------------
    cursors = [
        series_panel(ax_r1, movie, reduced, "rmsd_ca_tm", rmsd_ca, "CA RMSD (A)",
                     "TM CA RMSD vs the staged receptor  [02.11.00]",
                     "CA RMSD vs movie frame 0  [this movie]"),
        series_panel(ax_r2, movie, reduced, "lig_rmsd_pose", rmsd_lig, "ligand RMSD (A)",
                     "ligand RMSD vs the deposited pose  [02.11.00]",
                     "ligand RMSD vs movie frame 0  [this movie]"),
        series_panel(ax_r3, movie, reduced, "thickness", thickness, "P-P thickness (A)",
                     "bilayer thickness  [02.11.00]", "bilayer thickness  [this movie]"),
        series_panel(ax_r4, movie, reduced, "na_d250_dist", np.full(movie.n_frames, np.nan),
                     "Na+ - D2.50 (A)", "Na+ at the D2.50 site  [02.11.00]",
                     "Na+ - D2.50: needs the 02.11.00 reduction"),
        strip_cursor,
    ]
    cursors = [c for c in cursors if c is not None]

    aligned = movie.meta["movie"]["aligned"]
    header = fig.suptitle("", fontsize=11, x=0.035, ha="left", y=0.975)
    subhead = fig.text(0.035, 0.935, "", fontsize=7, ha="left", color="#555555")
    subhead.set_text(
        f"{movie.meta['atoms']['total']} atoms | "
        f"{movie.meta['movie']['n_frames']} frames every "
        f"{movie.meta['movie']['frame_spacing_ns']:.2f} ns | "
        f"residue numbering: {movie.meta['atoms']['numbering']} | "
        + (f"superposed on {movie.meta['movie']['align_selection']}" if aligned
           else "NOT superposed -- box drift, tilt and PBC artefacts are left in"))

    ca_ref = movie.xyz[0, movie.ca_ix, :]

    def update(i: int):
        pos = movie.xyz[i]
        ca = pos[movie.ca_ix]
        disp = np.linalg.norm(ca - ca_ref, axis=1)

        mem_lipid.set_offsets(pos[lipid_ix][:, [0, 2]])
        mem_phos.set_offsets(pos[phos_ix][:, [0, 2]])
        mem_water.set_offsets(pos[water_ix][:, [0, 2]])
        mem_ion.set_offsets(pos[ion_ix][:, [0, 2]])
        rec_seg = ca_segments(ca[is_rec], movie.chain_of_ca[is_rec], 0, 2)
        mem_rec.set_segments(rec_seg)
        # One colour per SEGMENT: the mean displacement of the two residues it joins.
        d_rec = disp[is_rec]
        same = movie.chain_of_ca[is_rec][:-1] == movie.chain_of_ca[is_rec][1:]
        mem_rec.set_array(0.5 * (d_rec[:-1][same] + d_rec[1:][same]))
        mem_oth.set_segments(ca_segments(ca[~is_rec], movie.chain_of_ca[~is_rec], 0, 2))
        mem_lig_pts.set_offsets(pos[lig_ix][:, [0, 2]] if has_ligand else np.empty((0, 2)))
        if lig_bonds.size:
            mem_lig.set_segments(np.stack([pos[lig_bonds[:, 0]][:, [0, 2]],
                                           pos[lig_bonds[:, 1]][:, [0, 2]]], axis=1))

        # Pocket view: a slab around the ligand, or the whole thing would overlap into mud.
        z0 = pos[centre_ix][:, 2].mean()
        def slab(ix):
            if not len(ix):
                return np.empty((0, 2))
            sel = ix[np.abs(pos[ix][:, 2] - z0) < POCKET_SLAB]
            return pos[sel][:, [0, 1]]
        pk_lipid.set_offsets(slab(lipid_ix))
        pk_water.set_offsets(slab(water_ix))
        pk_ion.set_offsets(slab(ion_ix))
        near = np.abs(ca[:, 2] - z0) < POCKET_SLAB
        pk_trace.set_segments(ca_segments(ca[near], movie.chain_of_ca[near], 0, 1))
        if contact_resids:
            pts = np.array([pos[movie.key[r]].mean(axis=0)[:2] for r in contact_resids])
            pk_key.set_offsets(pts)
            pk_key.set_facecolor([C["contact"] if contacts[i, j] else C["idle"]
                                  for j in range(len(contact_resids))])
            for t, p, r in zip(pk_labels, pts, contact_resids, strict=False):
                t.set_position((p[0], p[1] + 0.6))
                t.set_text(movie.key_labels[r])
        if lig_bonds.size:
            pk_lig.set_segments(np.stack([pos[lig_bonds[:, 0]][:, [0, 1]],
                                          pos[lig_bonds[:, 1]][:, [0, 1]]], axis=1))

        for cur in cursors:
            cur.set_xdata([movie.time_ns[i], movie.time_ns[i]])
        header.set_text(f"{movie.system} / {movie.replica}    "
                        f"t = {movie.time_ns[i]:7.1f} ns    "
                        f"frame {i + 1}/{movie.n_frames}")
        return ()

    return fig, update


def render(json_path: Path, out_dir: Path, fps: int, dpi: int,
           analysis_root: Path) -> Path:
    movie = Movie(json_path)
    reduced = reduced_series(movie.system, movie.replica, analysis_root)
    fig, update = build_figure(movie, reduced)
    fig.set_dpi(dpi)
    out = out_dir / (f"{PREFIX}_{movie.system}_{movie.replica}_dashboard_"
                     f"{date.today():%Y%m%d}.mp4")
    writer = pick_writer(fps)
    if writer is None:
        out = out.with_suffix(".gif")
        anim = animation.FuncAnimation(fig, update, frames=movie.n_frames, interval=1000 / fps)
        anim.save(str(out), writer=animation.PillowWriter(fps=fps), dpi=dpi)
    else:
        with writer.saving(fig, str(out), dpi):
            for i in range(movie.n_frames):
                update(i)
                writer.grab_frame()
    plt.close(fig)
    note = "" if reduced else "  (no 02.11.00 reduction found; traces come from the movie itself)"
    print(f"wrote {out}  [{movie.n_frames} frames, {movie.n_frames / fps:.1f} s]{note}")
    return out


def pick_writer(fps: int):
    """An ffmpeg writer, or None to fall back to an animated GIF."""
    if shutil.which("ffmpeg") is None:
        # Degrading rather than failing is deliberate -- a GIF is still watchable -- but it is
        # easy not to notice until the files are written, so the message names the fix rather
        # than the missing package.
        print("WARNING: ffmpeg is not on PATH; writing an animated GIF instead (several times "
              "larger, and not seekable). ffmpeg is declared in both environment files:\n"
              "    conda env update -f environment_zh853mor-local.yml   # local\n"
              "    conda env update -f environment_zh853mor-prep.yml    # cluster",
              file=sys.stderr)
        return None
    # -pix_fmt yuv420p so QuickTime and PowerPoint will play it; libx264 at CRF 20 keeps thin
    # lines (the CA trace) from turning to mush. The scale filter rounds the canvas down to even
    # dimensions: yuv420p cannot encode an odd width or height, and a non-default --dpi easily
    # produces one (e.g. 16.0 x 111 = 1776 x 999), which otherwise fails deep inside ffmpeg.
    return animation.FFMpegWriter(
        fps=fps, codec="libx264", bitrate=-1,
        extra_args=["-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
                    "-pix_fmt", "yuv420p", "-crf", "20", "-preset", "medium"])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--movies", type=Path, default=MOVIE_ROOT,
                    help="directory of exported movie trajectories")
    ap.add_argument("--analysis", type=Path, default=ANALYSIS_ROOT,
                    help="02.11.00 reduction directory, for the full-resolution traces")
    ap.add_argument("--system", help="only this system, e.g. apo_ASH")
    ap.add_argument("--replica", help="only this replica, e.g. prod_r1")
    ap.add_argument("--list", action="store_true",
                    help="print the numbered (system, replica) pairs and exit -- this is what "
                         "submit_render.sh sizes its job array from")
    ap.add_argument("--index", type=int,
                    help="render the Nth pair, 1-based -- what the SLURM array task passes")
    ap.add_argument("--fps", type=int, default=FPS)
    ap.add_argument("--dpi", type=int, default=DPI)
    ap.add_argument("--out", type=Path, default=paths.PRODUCT)
    args = ap.parse_args()

    found = sorted(args.movies.glob("*/*_movie.json"))
    if args.system:
        found = [p for p in found if p.parent.name == args.system]
    if args.replica:
        found = [p for p in found if p.name.startswith(f"{args.replica}_movie")]
    if not found:
        raise SystemExit(
            f"ERROR: no exported movies under {args.movies}. Run 01_export_movie_trajectory.py "
            "first (submit_export.sh on the cluster).")

    # One listing, used both to size the array and to resolve each task's index, so the mapping
    # cannot drift between the two -- the same arrangement 02.11.00 and the export step use.
    if args.list:
        for i, js in enumerate(found, start=1):
            meta = json.loads(js.read_text())
            mv = meta["movie"]
            print("\t".join([str(i), meta["build"], meta["replica"], str(mv["n_frames"]),
                             f"{mv['frame_spacing_ns']:.2f}", str(meta["atoms"]["total"])]))
        return 0
    if args.index is not None:
        if not 1 <= args.index <= len(found):
            raise SystemExit(f"ERROR: --index {args.index} outside 1-{len(found)}.")
        found = [found[args.index - 1]]

    paths.ensure_dir(args.out)
    failures = 0
    for js in found:
        try:
            render(js, args.out, args.fps, args.dpi, args.analysis)
        except Exception as exc:  # noqa: BLE001 -- one bad replica must not stop the panel
            print(f"FAILED {js.parent.name}/{js.stem}: {type(exc).__name__}: {exc}",
                  file=sys.stderr)
            failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

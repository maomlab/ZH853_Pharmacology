"""Time-lagged independent component analysis, and free-energy surfaces on the result.

Written out rather than taken from PyEMMA/deeptime for the reason given in D-9 for the
interaction fingerprint: the whole estimator is forty lines of linear algebra, the choices inside
it (which covariance is symmetrized, where the rank is truncated, how the sign is fixed) change
the answer, and they should be visible in this repository rather than in a dependency's defaults.
It also keeps the local analysis environment to the packages the rest of the project already uses.

**What tICA is for here.** PCA finds the directions of largest VARIANCE; tICA finds the directions
of slowest DECORRELATION (Molgedey & Schuster 1994; Perez-Hernandez 2013; Schwantes & Pande 2013).
For a receptor those are not the same thing: the largest-amplitude motion in a 500 ns window is
usually a floppy loop, while the interesting coordinate -- the ligand leaving a subsite, TM6
opening -- is slower and smaller. A landscape drawn on tICs therefore separates states that a
landscape drawn on PCs smears together.

**The trap it brings with it.** On a trajectory that has not converged, the slowest apparent
process IS the diffusive drift of an unequilibrated coordinate, so tIC1 reproduces the half-cosine
of free diffusion exactly as PC1 does (D-20). `zh853mor.convergence.cosine_content` applies to a
tIC projection unchanged, and the landscape stage reports it for the same reason.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Boltzmann constant in kcal/mol/K, and the production thermostat temperature (03_production.py).
KB_KCAL = 0.0019872041
TEMPERATURE_K = 310.0
KT_KCAL = KB_KCAL * TEMPERATURE_K          # 0.616 kcal/mol at 310 K


@dataclass(frozen=True)
class TICA:
    """A fitted tICA basis: what to subtract, what to project onto, and how slow each mode is."""

    mean: np.ndarray            # (n_features,) from the symmetrized ensemble, not from x_t alone
    components: np.ndarray      # (n_features, n_components), columns are the tICs
    eigenvalues: np.ndarray     # (n_components,) autocorrelation of each tIC at the fitted lag
    lag: int                    # frames
    rank: int                   # features kept after the whitening truncation
    n_samples: int              # time-lagged PAIRS used, not frames

    def transform(self, x: np.ndarray) -> np.ndarray:
        """Project features onto the tICs. `x` is (n_frames, n_features)."""
        x = np.asarray(x, dtype=float)
        if x.ndim != 2 or x.shape[1] != self.mean.size:
            raise ValueError(f"expected (n_frames, {self.mean.size}), got {x.shape}")
        return (x - self.mean) @ self.components

    def timescales(self, dt_ns: float = 1.0) -> np.ndarray:
        """Implied timescale of each tIC: -tau / ln(lambda), in units of `dt_ns` per frame.

        Only defined for 0 < lambda < 1, and the two ways out of that range mean OPPOSITE things,
        so they are reported differently rather than both collapsed to one sentinel:

        * **lambda >= 1** -- the mode has not decorrelated at all within the lag. Its timescale is
          longer than this trajectory can measure: `inf`.
        * **lambda <= 0** -- the mode has fully decorrelated and the residual sign is noise. Its
          timescale is SHORTER than the lag, not longer: `nan`. (Reporting `inf` here, as a naive
          guard does, turns the fastest modes in the model into the slowest.)
        """
        lam = np.asarray(self.eigenvalues, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            out = np.where((lam > 0) & (lam < 1), -self.lag * dt_ns / np.log(lam), np.inf)
        return np.where(lam <= 0, np.nan, out)


def _symmetrized_covariances(trajs: list[np.ndarray], lag: int
                             ) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Mean, instantaneous and time-lagged covariance, symmetrized over the pair (x_t, x_t+tau).

    Symmetrizing is what makes the estimator REVERSIBLE -- it enforces detailed balance, which is
    what licenses reading the eigenvalues as autocorrelations and the timescales as relaxation
    times. Without it C_tau is not symmetric, the eigenvalues can come out complex, and the
    "slowest mode" is partly a statement about the direction of time in a finite sample.

    Pairs are formed WITHIN a trajectory only. Pairing the end of one replica with the start of
    the next would manufacture a transition that never happened, and with three replicas that
    fiction lands squarely on the slowest timescale.
    """
    n = trajs[0].shape[1]
    c0 = np.zeros((n, n))
    ct = np.zeros((n, n))
    total = np.zeros(n)
    count = 0
    for x in trajs:
        if x.shape[0] <= lag:
            continue
        a, b = x[:-lag], x[lag:]
        c0 += a.T @ a + b.T @ b
        ct += a.T @ b + b.T @ a
        total += a.sum(axis=0) + b.sum(axis=0)
        count += a.shape[0]
    if count == 0:
        raise ValueError(f"no trajectory is longer than the lag of {lag} frames")
    mean = total / (2 * count)
    outer = np.outer(mean, mean)
    return mean, c0 / (2 * count) - outer, ct / (2 * count) - outer, count


def fit_tica(trajs: list[np.ndarray], lag: int, n_components: int = 2,
             epsilon: float = 1e-8, kinetic_map: bool = True) -> TICA:
    """Fit a tICA basis on one or more trajectories of features.

    `trajs` are (n_frames, n_features) arrays that must share their feature axis -- for a
    comparison across systems, fit on the POOLED list and project each system separately, or the
    axes of the two landscapes are different linear combinations and the pictures cannot be laid
    beside each other.

    `epsilon` truncates the whitening: features that are constant, duplicated, or linearly
    dependent (any pairwise-distance set is nearly so) give C0 near-singular eigenvalues, and
    dividing by their square roots turns rounding error into a spurious slowest mode. Directions
    below `epsilon` x the largest eigenvalue are dropped, and how many survived is reported as
    `rank` so the truncation is visible rather than implicit.

    With `kinetic_map` (Noe & Clementi 2015) each tIC is scaled by its eigenvalue, so Euclidean
    distance in the projection approximates kinetic distance and a histogram of it groups
    conformations that interconvert quickly. Turn it off to read the tICs as plain coordinates.
    """
    trajs = [np.asarray(t, dtype=float) for t in trajs]
    if not trajs:
        raise ValueError("no trajectories given")
    widths = {t.shape[1] for t in trajs}
    if len(widths) != 1:
        raise ValueError(f"trajectories have different feature counts: {sorted(widths)}")
    if lag < 1:
        raise ValueError("lag must be at least one frame")

    mean, c0, ct, count = _symmetrized_covariances(trajs, lag)

    # Whiten with C0, then diagonalize the lagged covariance in the whitened basis, where the
    # generalized problem C_tau v = lambda C0 v becomes an ordinary symmetric one.
    s, u = np.linalg.eigh(c0)
    keep = s > epsilon * max(s.max(), 0.0)
    if not keep.any():
        raise ValueError("the feature covariance is numerically zero; are all features constant?")
    whiten = u[:, keep] / np.sqrt(s[keep])
    k = whiten.T @ ct @ whiten
    k = 0.5 * (k + k.T)                      # kill the asymmetry left by floating-point error
    lam, vecs = np.linalg.eigh(k)
    order = np.argsort(lam)[::-1]            # slowest first
    lam, vecs = lam[order], vecs[:, order]

    n_components = min(n_components, vecs.shape[1])
    lam = lam[:n_components]
    comp = whiten @ vecs[:, :n_components]
    # The sign of an eigenvector is arbitrary, so fix it deterministically -- otherwise the same
    # data refitted flips an axis and two landscapes that should be identical are mirror images.
    flip = np.where(comp[np.argmax(np.abs(comp), axis=0), np.arange(n_components)] < 0, -1.0, 1.0)
    comp = comp * flip
    if kinetic_map:
        comp = comp * np.clip(lam, 0.0, None)
    return TICA(mean=mean, components=comp, eigenvalues=lam, lag=int(lag),
                rank=int(keep.sum()), n_samples=int(count))


def implied_timescale_scan(trajs: list[np.ndarray], lags: list[int], n_components: int = 4,
                           **kw) -> tuple[np.ndarray, np.ndarray]:
    """Implied timescales (in frames) against lag, the standard test that a lag is long enough.

    A timescale estimated at lag tau is only meaningful once it stops depending on tau: below
    that, the projection is still resolving processes faster than the lag and the estimate rises
    with it. Returns (lags, timescales[lag, component]); a fit that fails at a given lag -- every
    replica shorter than it -- contributes a row of NaN rather than dropping the lag silently.
    """
    out = np.full((len(lags), n_components), np.nan)
    ok = []
    for i, lag in enumerate(lags):
        try:
            model = fit_tica(trajs, lag=lag, n_components=n_components, **kw)
        except (ValueError, np.linalg.LinAlgError):
            continue
        ts = model.timescales()
        out[i, :ts.size] = ts
        ok.append(lag)
    return np.asarray(lags, dtype=float), out


def common_extent(arrays: list[np.ndarray], pad: float = 0.05) -> tuple[float, float]:
    """(lo, hi) covering every array, padded -- so every system is binned on the SAME grid."""
    finite = np.concatenate([np.asarray(a, dtype=float).ravel() for a in arrays])
    finite = finite[np.isfinite(finite)]
    if not finite.size:
        return (0.0, 1.0)
    lo, hi = float(finite.min()), float(finite.max())
    span = hi - lo or 1.0
    return lo - pad * span, hi + pad * span


def free_energy(x: np.ndarray, y: np.ndarray, xedges: np.ndarray, yedges: np.ndarray,
                kt: float = KT_KCAL, min_count: int = 5) -> tuple[np.ndarray, np.ndarray]:
    """-kT ln P on a fixed grid, zeroed at its own minimum. Returns (free energy, counts).

    This is a SAMPLING density turned into energy units, not a converged free energy (D-5): the
    depth of a basin the trajectory visited twice is a statement about those two visits. Bins with
    fewer than `min_count` frames are returned as NaN rather than as a very high barrier, because
    -kT ln(1/N) for a single stray frame draws a confident-looking wall out of one sample.

    Zeroing at the minimum makes the colour scale read as "cost relative to the most populated
    state", which is comparable between systems only because the caller bins them identically.
    """
    counts, _, _ = np.histogram2d(np.asarray(x, dtype=float), np.asarray(y, dtype=float),
                                  bins=[xedges, yedges])
    total = counts.sum()
    if total <= 0:
        return np.full(counts.shape, np.nan), counts
    p = counts / total
    with np.errstate(divide="ignore", invalid="ignore"):
        f = -kt * np.log(p)
    f[counts < min_count] = np.nan
    if np.isfinite(f).any():
        f -= np.nanmin(f)
    return f, counts

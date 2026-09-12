"""Convergence and uncertainty statistics for MD observables.

The Phase-4 gate in ``docs/PLAN.md`` asks for "a documented convergence assessment (not
RMSD-plateau alone)". A plateau is necessary and nowhere near sufficient: a correlated series
that has stopped drifting still carries far fewer independent samples than it has frames, and
averaging over them produces error bars that are too small by exactly the factor this module
computes.

Four complementary diagnostics, none of which can be read off a time series by eye:

``statistical_inefficiency``/``detect_equilibration``
    How many frames one independent sample costs (g), where the production part of the series
    starts, and therefore the effective sample size. Chodera's marginal-effective-sample-size
    rule (J. Chem. Theory Comput. 2016, 12, 1799) chooses the discard point by maximising the
    number of uncorrelated samples that survive it, rather than by eye.
``block_sem_curve``
    The textbook blocking check: the standard error estimated from block means rises with block
    size and plateaus once blocks are longer than the correlation time. A curve that never
    plateaus means the run is shorter than its own correlation time and the error bar is a lower
    bound, whatever ``statistical_inefficiency`` reports.
``gelman_rubin``
    Do the independent replicas agree? Within-replica variance alone cannot see a system stuck
    in one basin; comparing it against the between-replica variance can. R-hat near 1 means the
    replicas are sampling one distribution.
``cosine_content``
    Hess's test (Phys. Rev. E 2002, 65, 031910) for the specific failure where a principal
    component of a random walk imitates a slow collective motion. The projection of free
    diffusion onto its own first PC is a half-cosine, so a cosine content near 1 says the
    "motion" is diffusive noise and near 0 that it is real sampling.

Pure NumPy, no MDAnalysis: this is the half of the analysis that is unit-testable against series
with known answers, and ``tests/test_convergence.py`` does exactly that.
"""

from __future__ import annotations

import numpy as np


def _as_1d(x: np.ndarray | list[float]) -> np.ndarray:
    a = np.asarray(x, dtype=float).ravel()
    if a.size == 0:
        raise ValueError("empty series")
    return a


def autocorrelation(x: np.ndarray | list[float]) -> np.ndarray:
    """Normalised autocorrelation function C(t), C(0) = 1, via FFT.

    Unbiased in the sense that each lag is divided by the number of pairs contributing to it,
    so long lags are noisy rather than damped towards zero -- which is why every consumer here
    truncates at the first non-positive value instead of summing the whole thing.
    """
    a = _as_1d(x)
    n = a.size
    a = a - a.mean()
    var = float(np.dot(a, a))
    if var <= 0:  # a constant series is perfectly correlated with itself and nothing else
        return np.ones(1)
    size = 1 << int(np.ceil(np.log2(2 * n)))
    f = np.fft.rfft(a, size)
    acf = np.fft.irfft(f * np.conjugate(f), size)[:n]
    return acf / (var * (1.0 - np.arange(n) / n))


def statistical_inefficiency(x: np.ndarray | list[float]) -> float:
    """g = 1 + 2 tau: frames per independent sample (>= 1).

    N/g is the effective number of independent samples, so the standard error of the mean is
    std/sqrt(N/g) and NOT std/sqrt(N). For a 500 ns trajectory sampled every 100 ps, g is
    routinely 10-100 for pocket observables: ignoring it overstates the precision tenfold.
    """
    a = _as_1d(x)
    n = a.size
    if n < 3 or np.allclose(a, a[0]):
        return 1.0
    c = autocorrelation(a)
    g = 1.0
    for t in range(1, n):
        if c[t] <= 0:  # first zero crossing: beyond it the ACF is noise, not signal
            break
        g += 2.0 * (1.0 - t / n) * c[t]
    return float(max(g, 1.0))


def detect_equilibration(x: np.ndarray | list[float], n_candidates: int = 100) -> tuple[int, float, float]:
    """(t0, g, n_eff): where production starts, by maximum effective sample size.

    Scans candidate discard points and keeps the one leaving the most uncorrelated samples
    behind it (Chodera 2016). Discarding too little leaves the relaxation in the average;
    discarding too much throws away the sampling -- this trades the two off explicitly instead
    of applying a fixed "first 10%" rule that is wrong in both directions on different systems.
    """
    a = _as_1d(x)
    n = a.size
    if n < 10:
        return 0, statistical_inefficiency(a), float(n)
    # Candidates over the first 80%: a t0 past that leaves too little series to characterise.
    stop = max(int(0.8 * n), 1)
    step = max(stop // n_candidates, 1)
    best = (0, 1.0, -np.inf)
    for t0 in range(0, stop, step):
        tail = a[t0:]
        if tail.size < 3:
            break
        g = statistical_inefficiency(tail)
        n_eff = tail.size / g
        if n_eff > best[2]:
            best = (t0, g, n_eff)
    return int(best[0]), float(best[1]), float(best[2])


def mean_sem(x: np.ndarray | list[float], g: float | None = None) -> tuple[float, float]:
    """Mean and correlation-corrected standard error, std/sqrt(N/g)."""
    a = _as_1d(x)
    if a.size < 2:
        return float(a[0]), float("nan")
    if g is None:
        g = statistical_inefficiency(a)
    n_eff = max(a.size / max(g, 1.0), 1.0)
    return float(a.mean()), float(a.std(ddof=1) / np.sqrt(n_eff))


def block_sem_curve(x: np.ndarray | list[float], max_blocks: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """(block sizes, SEM estimated from block means) for a blocking analysis.

    Block sizes double from 1 up to N/`max_blocks` blocks (at least 2 blocks per point, or the
    SEM of the block means is undefined). The SEM plateau value is the honest error bar; the
    presence of a plateau is the convergence statement.
    """
    a = _as_1d(x)
    n = a.size
    sizes, sems = [], []
    size = 1
    while n // size >= max(max_blocks, 2):
        nb = n // size
        means = a[: nb * size].reshape(nb, size).mean(axis=1)
        sizes.append(size)
        sems.append(float(means.std(ddof=1) / np.sqrt(nb)))
        size *= 2
    return np.array(sizes, dtype=int), np.array(sems, dtype=float)


def gelman_rubin(chains: list[np.ndarray] | np.ndarray) -> float:
    """Potential scale reduction R-hat across independent replicas.

    R-hat ~ 1.0 means the replicas cannot be told apart from one distribution; > 1.1 is the
    conventional "not converged" flag. With the 3 replicas of this study R-hat is a coarse
    instrument -- it is reported as a flag, not as a precise number.
    """
    arrs = [_as_1d(c) for c in chains]
    m = len(arrs)
    if m < 2:
        return float("nan")
    n = min(a.size for a in arrs)
    if n < 2:
        return float("nan")
    a = np.array([c[-n:] for c in arrs])  # common length, taking the LATEST n of each
    means = a.mean(axis=1)
    within = float(a.var(axis=1, ddof=1).mean())
    between = float(n * means.var(ddof=1))
    if within <= 0:
        return float("nan")
    var_hat = (n - 1) / n * within + between / n
    return float(np.sqrt(var_hat / within))


def cosine_content(projection: np.ndarray | list[float], index: int = 1) -> float:
    """Hess cosine content of a PCA projection: 1 = free diffusion, 0 = real sampling.

    A high value on the first principal components means PCA has fitted the half-cosine shape
    of a random walk, so the apparent large-amplitude "collective motion" carries no
    conformational information.
    """
    p = _as_1d(projection)
    n = p.size
    denom = float(np.dot(p, p))
    if denom <= 0:
        return 0.0
    t = np.arange(n)
    num = float(np.dot(np.cos(index * np.pi * t / n), p))
    return float(2.0 / n * num**2 / denom)


def linear_drift(y: np.ndarray | list[float], x: np.ndarray | list[float] | None = None) -> tuple[float, float]:
    """(slope, standard error of slope) from an ordinary least-squares fit.

    The plateau test with a number attached: a drift whose magnitude is within ~2 standard
    errors of zero is not resolvable, which is a statement a "looks flat" verdict cannot make.
    The SE is NOT corrected for autocorrelation, so it is optimistic; treat it as a screen.
    """
    a = _as_1d(y)
    t = np.arange(a.size, dtype=float) if x is None else _as_1d(x)
    if a.size < 3:
        return float("nan"), float("nan")
    slope, intercept = np.polyfit(t, a, 1)
    resid = a - (slope * t + intercept)
    ss = float(np.dot(resid, resid)) / (a.size - 2)
    sxx = float(np.dot(t - t.mean(), t - t.mean()))
    return float(slope), float(np.sqrt(ss / sxx)) if sxx > 0 else float("nan")

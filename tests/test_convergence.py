"""Tests for zh853mor.convergence against series whose answers are known analytically.

The point of these statistics is that a time series cannot be judged by eye, so they are
validated against processes with a closed-form answer rather than against a trajectory.
"""

from __future__ import annotations

import numpy as np
import pytest

from zh853mor import convergence as cv


def ar1(n: int, phi: float, seed: int = 0) -> np.ndarray:
    """AR(1) process: known statistical inefficiency g = (1 + phi) / (1 - phi)."""
    rng = np.random.default_rng(seed)
    noise = rng.normal(size=n)
    x = np.empty(n)
    x[0] = noise[0]
    for i in range(1, n):
        x[i] = phi * x[i - 1] + noise[i]
    return x


def test_uncorrelated_series_has_g_near_one():
    x = np.random.default_rng(1).normal(size=20000)
    assert cv.statistical_inefficiency(x) == pytest.approx(1.0, abs=0.5)


def test_ar1_recovers_analytic_g():
    """g = (1+phi)/(1-phi) = 9 for phi = 0.8 -- i.e. 9 frames per independent sample."""
    g = cv.statistical_inefficiency(ar1(200000, 0.8))
    assert g == pytest.approx(9.0, rel=0.25)


def test_constant_and_short_series_are_safe():
    assert cv.statistical_inefficiency(np.ones(100)) == 1.0
    assert cv.statistical_inefficiency([1.0, 2.0]) == 1.0


def test_correlated_sem_exceeds_naive_sem():
    x = ar1(20000, 0.9)
    _, sem = cv.mean_sem(x)
    naive = x.std(ddof=1) / np.sqrt(x.size)
    assert sem > 3 * naive  # sqrt(g) ~ 4.4; ignoring correlation would overstate precision


def test_detect_equilibration_discards_a_transient():
    """An exponential relaxation followed by a stationary tail: t0 lands inside the transient."""
    n, tau = 4000, 200.0
    t = np.arange(n)
    x = 10.0 * np.exp(-t / tau) + np.random.default_rng(3).normal(size=n)
    t0, g, n_eff = cv.detect_equilibration(x)
    assert 0 < t0 < 2000
    assert n_eff <= n - t0
    assert g >= 1.0
    # The mean after discarding is close to the stationary value; before, it is not.
    assert abs(float(x[t0:].mean())) < abs(float(x.mean()))


def test_block_sem_curve_plateaus_at_the_correlated_sem():
    x = ar1(20000, 0.9)
    sizes, sems = cv.block_sem_curve(x)
    assert sizes[0] == 1
    assert sems[-1] > sems[0]                      # blocking reveals the hidden correlation
    _, sem = cv.mean_sem(x)
    assert sems[-3:].mean() == pytest.approx(sem, rel=0.6)


def test_gelman_rubin_flags_disagreeing_replicas():
    rng = np.random.default_rng(5)
    same = [rng.normal(size=2000) for _ in range(3)]
    assert cv.gelman_rubin(same) == pytest.approx(1.0, abs=0.05)
    shifted = [rng.normal(size=2000) + offset for offset in (0.0, 1.0, 2.0)]
    assert cv.gelman_rubin(shifted) > 1.1
    assert np.isnan(cv.gelman_rubin([np.zeros(10)]))


def test_cosine_content_separates_diffusion_from_sampling():
    n = 5000
    t = np.arange(n)
    assert cv.cosine_content(np.cos(np.pi * t / n)) == pytest.approx(1.0, abs=0.01)
    assert cv.cosine_content(np.random.default_rng(7).normal(size=n)) < 0.05
    # The first PC of a random walk is the case this exists to catch.
    walk = np.cumsum(np.random.default_rng(11).normal(size=n))
    assert cv.cosine_content(walk - walk.mean()) > 0.5


def test_linear_drift_recovers_a_known_slope():
    n = 1000
    t = np.arange(n, dtype=float)
    y = 0.01 * t + np.random.default_rng(13).normal(scale=0.1, size=n)
    slope, se = cv.linear_drift(y)
    assert slope == pytest.approx(0.01, rel=0.05)
    assert se < 0.001
    # Calibration of the error bar on the null: over many pure-noise series, the slope should
    # sit within 2 SE of zero about 95% of the time. Asserting it for one fixed draw would be a
    # test of the seed, not of the fit.
    rng = np.random.default_rng(17)
    within = 0
    for _ in range(100):
        s2, se2 = cv.linear_drift(rng.normal(size=n))
        within += abs(s2) < 2 * se2
    assert within >= 88

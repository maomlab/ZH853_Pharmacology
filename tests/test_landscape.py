"""Tests for zh853mor.landscape on processes whose slow mode is known by construction.

tICA is checked against an Ornstein-Uhlenbeck process, because that is the one case where the
right answer is available in closed form: an OU coordinate with relaxation time tau has
autocorrelation exp(-lag/tau), so the implied timescale the estimator reports must come back as
tau. Mixing two such coordinates by a rotation then asks the real question -- can it find the slow
direction when no single input feature is it, which is the entire reason for using tICA over
watching individual distances.
"""

from __future__ import annotations

import numpy as np
import pytest

from zh853mor import convergence as cv
from zh853mor import landscape as ls

TAU_SLOW, TAU_FAST = 200.0, 5.0        # frames
N_FRAMES = 60_000                      # long enough that a 200-frame process is sampled ~300x


def _ou(tau: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """Stationary OU path with relaxation time `tau` and unit variance."""
    a = np.exp(-1.0 / tau)
    noise = rng.normal(0.0, np.sqrt(1.0 - a * a), n)
    x = np.empty(n)
    x[0] = rng.normal()
    for i in range(1, n):
        x[i] = a * x[i - 1] + noise[i]
    return x


def _mixed(rng, n=N_FRAMES, angle=0.7, amp_fast=3.0):
    """Slow and fast OU coordinates, rotated together and with the FAST one made larger.

    The fast coordinate carries the bigger amplitude on purpose: PCA would return it first, so a
    test that passes here is testing tICA rather than a direction any variance-based method would
    have found anyway.
    """
    slow, fast = _ou(TAU_SLOW, n, rng), amp_fast * _ou(TAU_FAST, n, rng)
    c, s = np.cos(angle), np.sin(angle)
    rot = np.array([[c, -s], [s, c]])
    return np.stack([slow, fast], axis=1) @ rot.T, rot[:, 0]


def test_tica_finds_the_slow_direction_that_pca_would_miss():
    rng = np.random.default_rng(0)
    x, slow_dir = _mixed(rng)

    model = ls.fit_tica([x], lag=50, n_components=2, kinetic_map=False)
    tic1 = model.components[:, 0] / np.linalg.norm(model.components[:, 0])
    assert abs(float(tic1 @ slow_dir)) > 0.97

    # And confirm the premise: the leading PCA direction is the loud fast one, not this.
    cov = np.cov(x.T)
    pc1 = np.linalg.eigh(cov)[1][:, -1]
    assert abs(float(pc1 @ slow_dir)) < 0.5


def test_implied_timescale_recovers_the_generating_relaxation_time():
    rng = np.random.default_rng(1)
    x, _ = _mixed(rng)
    model = ls.fit_tica([x], lag=50, n_components=2, kinetic_map=False)
    t1, t2 = model.timescales()
    assert TAU_SLOW * 0.75 < t1 < TAU_SLOW * 1.35
    assert t2 < TAU_SLOW / 5          # the fast mode is resolved as clearly faster


def test_timescales_plateau_once_the_lag_is_long_enough():
    """The standard justification for a lag choice: the estimate stops depending on it."""
    rng = np.random.default_rng(2)
    x, _ = _mixed(rng)
    lags, ts = ls.implied_timescale_scan([x], lags=[2, 10, 50, 100, 200], n_components=2,
                                         kinetic_map=False)
    t1 = ts[:, 0]
    assert np.isfinite(t1).all()
    assert t1[0] < t1[2]                                        # too-short lag underestimates
    plateau = t1[2:]
    assert plateau.max() / plateau.min() < 1.6                  # and then it settles


def test_pairs_never_cross_a_trajectory_boundary():
    """Replicas must not be glued into one: the joint is a transition that never happened.

    Checked structurally -- the pair count is sum(T_i - lag), not (sum T_i) - lag -- and then
    behaviourally, on independent realizations of the SAME slow process. Pairs that straddle the
    joint relate two unrelated samples, so they drag the estimated autocorrelation down and the
    timescale with it. Short replicas and a long lag are used deliberately: that is the regime
    where the straddling pairs are a large enough fraction to see, and it is also the regime this
    project is actually in (3 x 500 ns replicas of a process worth tens of ns).
    """
    rng = np.random.default_rng(3)
    n, lag = 400, 60
    reps = [_ou(TAU_SLOW, n, rng).reshape(-1, 1) for _ in range(6)]

    separate = ls.fit_tica(reps, lag=lag, n_components=1, kinetic_map=False)
    glued = ls.fit_tica([np.vstack(reps)], lag=lag, n_components=1, kinetic_map=False)

    assert separate.n_samples == sum(len(r) - lag for r in reps)
    assert glued.n_samples == len(np.vstack(reps)) - lag
    assert glued.timescales()[0] < separate.timescales()[0]


def test_projection_of_a_drifting_coordinate_is_caught_by_cosine_content():
    """The D-20 trap, in tICA form: on unconverged data the slowest mode IS the drift."""
    rng = np.random.default_rng(4)
    n = 5000
    drift = np.linspace(0, 12, n) + 0.3 * rng.normal(size=n)     # never equilibrates
    noise = rng.normal(size=n)
    x = np.stack([drift, noise], axis=1)
    model = ls.fit_tica([x], lag=25, n_components=1, kinetic_map=False)
    assert cv.cosine_content(model.transform(x)[:, 0], index=1) > 0.5


def test_sign_is_fixed_so_refitting_cannot_mirror_the_landscape():
    rng = np.random.default_rng(5)
    x, _ = _mixed(rng, n=20_000)
    a = ls.fit_tica([x], lag=50, n_components=2)
    b = ls.fit_tica([x], lag=50, n_components=2)
    assert np.allclose(a.components, b.components)
    biggest = np.abs(a.components).argmax(axis=0)
    assert (a.components[biggest, np.arange(a.components.shape[1])] > 0).all()


def test_rank_truncation_survives_linearly_dependent_features():
    """Pairwise-distance features are always nearly dependent; an exact duplicate is the limit."""
    rng = np.random.default_rng(6)
    x, _ = _mixed(rng, n=20_000)
    padded = np.column_stack([x, x[:, 0], np.full(len(x), 2.5)])   # a copy and a constant
    model = ls.fit_tica([padded], lag=50, n_components=2, kinetic_map=False)
    assert model.rank == 2                                         # not 4
    assert TAU_SLOW * 0.75 < model.timescales()[0] < TAU_SLOW * 1.35


def test_a_fully_decorrelated_mode_is_nan_not_inf():
    """A negative eigenvalue means FASTER than the lag; reporting inf would invert its meaning."""
    rng = np.random.default_rng(12)
    x, _ = _mixed(rng, n=20_000)
    # lag 50 vs TAU_FAST 5: the fast mode is long gone, so its eigenvalue is zero plus noise.
    model = ls.fit_tica([x], lag=50, n_components=2, kinetic_map=False)
    ts = model.timescales()
    assert np.isfinite(ts[0])
    if model.eigenvalues[1] <= 0:
        assert np.isnan(ts[1])
    else:
        assert ts[1] < TAU_SLOW / 5


def test_kinetic_map_scales_each_component_by_its_eigenvalue():
    rng = np.random.default_rng(7)
    x, _ = _mixed(rng, n=20_000)
    plain = ls.fit_tica([x], lag=50, n_components=2, kinetic_map=False)
    scaled = ls.fit_tica([x], lag=50, n_components=2, kinetic_map=True)
    assert np.allclose(scaled.components, plain.components * np.clip(plain.eigenvalues, 0, None))


def test_a_lag_longer_than_every_trajectory_is_an_error_not_an_empty_answer():
    rng = np.random.default_rng(8)
    short = [rng.normal(size=(30, 4)), rng.normal(size=(20, 4))]
    with pytest.raises(ValueError, match="longer than the lag"):
        ls.fit_tica(short, lag=100)


def test_mismatched_feature_counts_are_refused():
    rng = np.random.default_rng(9)
    with pytest.raises(ValueError, match="different feature counts"):
        ls.fit_tica([rng.normal(size=(100, 4)), rng.normal(size=(100, 5))], lag=5)


def test_free_energy_is_zeroed_at_its_minimum_and_masks_thin_bins():
    rng = np.random.default_rng(10)
    x = rng.normal(size=20_000)
    y = rng.normal(size=20_000)
    xe = np.linspace(-4, 4, 21)
    ye = np.linspace(-4, 4, 21)
    f, counts = ls.free_energy(x, y, xe, ye, min_count=5)

    assert np.nanmin(f) == pytest.approx(0.0)
    assert np.isnan(f[counts < 5]).all()
    # A standard normal is a harmonic well: F(1 sigma) - F(0) = kT/2 per dimension.
    centre = f[9:11, 9:11]
    assert np.nanmean(centre) < 0.2
    assert np.isfinite(f[counts >= 5]).all()


def test_free_energy_bins_are_the_callers_so_systems_are_comparable():
    """Two different samples binned on one grid must produce alignable arrays."""
    rng = np.random.default_rng(11)
    a = rng.normal(size=(5000, 2))
    b = rng.normal(loc=1.5, size=(5000, 2))
    xe = np.linspace(*ls.common_extent([a[:, 0], b[:, 0]]), 25)
    ye = np.linspace(*ls.common_extent([a[:, 1], b[:, 1]]), 25)
    fa, _ = ls.free_energy(a[:, 0], a[:, 1], xe, ye)
    fb, _ = ls.free_energy(b[:, 0], b[:, 1], xe, ye)
    assert fa.shape == fb.shape
    diff = fa - fb
    assert np.isfinite(diff).any()      # they overlap, so a difference map is defined somewhere


# --- the 02.13.00 feature resolver ---------------------------------------------------------------
# `--space pair` lets a caller name a coordinate as a string. Every one of those strings indexes
# into an array by position, so a silently wrong lookup produces a plausible landscape of the
# wrong coordinate -- the failure mode with no symptom.

import importlib.util  # noqa: E402
import sys  # noqa: E402

from zh853mor import md, paths  # noqa: E402

_STEP = paths.SRC / "02.13.00_map_conformational_landscapes" / "01_fit_landscape.py"


def _fit_module():
    spec = importlib.util.spec_from_file_location("fit_landscape", _STEP)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def reduced(tmp_path):
    """A reduced replica shaped exactly as 01_reduce_trajectory.py writes one."""
    import json
    n, t0 = 50, 10
    d = tmp_path / "ZH853_ASP"
    d.mkdir()
    rng = np.random.default_rng(0)
    np.savez_compressed(
        d / "prod_r1.npz",
        time_ns=np.arange(n, dtype=np.float32) * 0.1,
        na_d250_dist=np.linspace(3.0, 5.0, n).astype(np.float32),
        activation=np.column_stack([np.linspace(8, 12, n), np.linspace(11, 9, n)]),
        activation_labels=np.array(["tm3_tm6_ca", "npxxy_tm3_ca"]),
        resids=np.array([149, 299, 328]),
        min_dist=rng.uniform(3, 9, (n, 3)).astype(np.float16),
        key_pair_ids=np.array([[116, 167], [116, 281], [167, 281]]),
        key_pairs=rng.uniform(6, 20, (n, 3)).astype(np.float16),
    )
    (d / "prod_r1.json").write_text(json.dumps(
        {"build": "ZH853_ASP", "replica": "prod_r1",
         "equilibration": {"t0_frames": t0, "t0_ns": 1.0}}))
    return md.load_replicas(tmp_path)[0], t0, n


def test_named_features_resolve_to_the_right_column(reduced):
    fit = _fit_module()
    result, t0, n = reduced
    arrays = result.arrays()

    plain = fit.named_feature(result, "na_d250_dist")
    assert plain.size == n - t0
    assert np.allclose(plain, arrays["na_d250_dist"][t0:])

    # The SECOND ruler, by name -- the case a positional lookup gets silently wrong.
    npxxy = fit.named_feature(result, "activation:npxxy_tm3_ca")
    assert np.allclose(npxxy, arrays["activation"][t0:, 1])

    assert np.allclose(fit.named_feature(result, "lig:299"),
                       np.asarray(arrays["min_dist"], dtype=float)[t0:, 1])
    assert np.allclose(fit.named_feature(result, "ca:116-281"),
                       np.asarray(arrays["key_pairs"], dtype=float)[t0:, 1])
    # A CA pair is unordered: naming it backwards must find the same column, not raise.
    assert np.allclose(fit.named_feature(result, "ca:281-116"),
                       fit.named_feature(result, "ca:116-281"))


def test_named_features_start_after_equilibration(reduced):
    """Fitting on the relaxation would make tIC1 the equilibration transient, every time."""
    fit = _fit_module()
    result, t0, n = reduced
    for name in ("na_d250_dist", "activation:tm3_tm6_ca", "lig:149", "ca:116-167"):
        assert fit.named_feature(result, name).size == n - t0


@pytest.mark.parametrize("name,match", [
    ("activation:nope", "not among"),
    ("lig:999", "not in the construct"),
    ("ca:116-999", "no CA-CA pair"),
    ("bogus:1", "unknown feature prefix"),
])
def test_unknown_features_are_named_not_guessed(reduced, name, match):
    fit = _fit_module()
    result, _, _ = reduced
    with pytest.raises(KeyError, match=match):
        fit.named_feature(result, name)


def test_between_system_fraction_separates_a_label_from_a_coordinate():
    """The number that says whether a tIC is a kinetic mode or just 'which simulation is this'."""
    fit = _fit_module()
    rng = np.random.default_rng(1)
    # Component 0 separates the systems completely; component 1 is the same in both.
    a = np.column_stack([rng.normal(-5, 0.3, 500), rng.normal(0, 1.0, 500)])
    b = np.column_stack([rng.normal(+5, 0.3, 500), rng.normal(0, 1.0, 500)])
    frac = fit.between_system_fraction({"A": [("r1", a, None, None)],
                                        "B": [("r1", b, None, None)]})
    assert frac[0] > 0.9
    assert frac[1] < 0.1


def test_between_system_fraction_is_zero_for_a_single_system():
    fit = _fit_module()
    rng = np.random.default_rng(2)
    one = {"A": [("r1", rng.normal(size=(200, 2)), None, None)]}
    assert np.allclose(fit.between_system_fraction(one), 0.0)

"""
Test suite for the CSD inference engine.

These tests double as the scientific validation: they assert that the tool
(1) recovers known OU parameters, (2) stays quiet on stationary noise,
(3) detects an approaching fold, and (4) rejects a non-OU series through its
goodness-of-fit gate. If these pass, the engine behaves as documented.

Run with:  pytest -q
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import csd_inference as ci


def _series(x, dt=1.0, name="s"):
    return ci.Series(np.arange(len(x), dtype=float), x, dt, name)


def test_ou_parameter_recovery():
    """Estimated k and sigma should bracket the true values within CI."""
    true_k, true_sigma = 0.2, 1.0
    x = ci.simulate_ou(2000, 1.0, k=true_k, sigma=true_sigma,
                       rng=np.random.default_rng(0))
    fit = ci.estimate_ou(_series(x), n_boot=200, rng=np.random.default_rng(0))
    assert fit.k_ci[0] <= true_k <= fit.k_ci[1], (true_k, fit.k_ci)
    assert fit.sigma_ci[0] <= true_sigma <= fit.sigma_ci[1], (true_sigma, fit.sigma_ci)


def test_stationary_control_not_significant():
    """A constant-k process must not produce a significant CSD trend."""
    x = ci.simulate_ou(800, 1.0, k=0.2, sigma=1.0,
                       rng=np.random.default_rng(2002))
    tt = ci.surrogate_significance(_series(x), n_surrogates=300,
                                   rng=np.random.default_rng(2))
    assert tt.p_ar1 > 0.05, tt.p_ar1


def test_approaching_fold_is_detected():
    """A decaying restoring rate must yield a significant, positive trend."""
    x = ci.simulate_approaching_fold(800, 1.0, 0.4, 0.02, 1.0,
                                     np.random.default_rng(1))
    tt = ci.surrogate_significance(_series(x), n_surrogates=300,
                                   rng=np.random.default_rng(1))
    assert tt.tau_ar1 > 0.5
    assert tt.p_ar1 < 0.05, tt.p_ar1


def test_goodness_of_fit_rejects_non_ou():
    """A strongly periodic (non-OU) series must fail the adequacy gate."""
    t = np.arange(800.0)
    osc = ci.simulate_ou(800, 1.0, k=0.3, sigma=0.5,
                         rng=np.random.default_rng(4)) + 3 * np.sin(2 * np.pi * t / 12)
    gof = ci.ou_goodness_of_fit(ci.Series(t, osc, 1.0))
    assert not gof.adequate
    assert gof.ljung_box_p < 0.05


def test_goodness_of_fit_accepts_ou():
    """A clean OU series must pass the adequacy gate."""
    x = ci.simulate_ou(800, 1.0, k=0.25, sigma=1.0,
                       rng=np.random.default_rng(7))
    gof = ci.ou_goodness_of_fit(_series(x))
    assert gof.adequate


def test_rolling_indicators_match_naive():
    """Vectorised rolling AR(1) must equal the direct per-window computation."""
    x = ci.simulate_ou(300, 1.0, k=0.2, sigma=1.0, rng=np.random.default_rng(3))
    ser = _series(x)
    r = ci.rolling_indicators(ser, 0.5)
    resid = ci.detrend(x, r.bandwidth)
    W = r.window
    for i in (0, 40, 90):
        naive = ci._lag1_autocorr(resid[i:i + W])
        assert abs(r.ar1[i] - naive) < 1e-9, (i, r.ar1[i], naive)


def test_permutation_null_weaker_than_fourier():
    """The permutation null destroys all autocorrelation (white noise), so it
    is weaker and should give a p-value <= the Fourier null on a correlated
    series with a trend. Both must be valid probabilities."""
    x = ci.simulate_approaching_fold(600, 1.0, 0.4, 0.03, 1.0,
                                     np.random.default_rng(1))
    ser = _series(x)
    f = ci.surrogate_significance(ser, n_surrogates=200, null_method="fourier",
                                  rng=np.random.default_rng(0))
    p = ci.surrogate_significance(ser, n_surrogates=200,
                                  null_method="permutation",
                                  rng=np.random.default_rng(0))
    assert 0 < f.p_ar1 <= 1 and 0 < p.p_ar1 <= 1
    assert p.p_ar1 <= f.p_ar1 + 1e-9        # permutation no stricter
    assert f.tau_ar1 == p.tau_ar1           # same observed statistic


def test_reproducibility():
    """Same seed => identical surrogate p-value."""
    x = ci.simulate_ou(400, 1.0, k=0.2, sigma=1.0, rng=np.random.default_rng(5))
    ser = _series(x)
    a = ci.surrogate_significance(ser, n_surrogates=200, rng=np.random.default_rng(9))
    b = ci.surrogate_significance(ser, n_surrogates=200, rng=np.random.default_rng(9))
    assert a.p_ar1 == b.p_ar1


# ---------------------------------------------------------------------------
#  Regression tests for real-instrument-data hardening (the ICU SpO2 case):
#  ambiguous timestamps must be refused, and sensor dropouts must be
#  cleanable so a buried CSD signal is recovered.
# ---------------------------------------------------------------------------
def _write_csv(path, t, x, cols=("time", "value")):
    import csv
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(list(cols))
        for ti, xi in zip(t, x):
            w.writerow([ti, float(xi)])


def test_loader_refuses_ambiguous_dt(tmp_path):
    """Rounded/duplicate timestamps (median diff 0) must raise, not silently
    default to dt=1.0. This is the failure the ICU run hit."""
    x = ci.simulate_ou(600, 1.0, k=0.2, sigma=1.0, rng=np.random.default_rng(0))
    # each timestamp repeated 3x => majority of consecutive diffs are 0 => median 0
    t = np.floor(np.arange(600) / 3.0)
    p = tmp_path / "rounded.csv"
    _write_csv(p, t, x)
    try:
        ci.load_series(str(p))
        assert False, "expected ValueError on ambiguous dt"
    except ValueError as e:
        assert "sampling interval" in str(e)


def test_dt_override_accepted(tmp_path):
    """An explicit dt must bypass timestamp inference entirely."""
    x = ci.simulate_ou(600, 1.0, k=0.2, sigma=1.0, rng=np.random.default_rng(0))
    t = np.floor(np.arange(600) / 3.0)
    p = tmp_path / "rounded.csv"
    _write_csv(p, t, x)
    ser = ci.load_series(str(p), dt_override=1.0)
    assert ser.dt == 1.0
    assert ser.n == 600


def test_dropout_cleaning_recovers_signal():
    """A genuine approaching-fold signal, contaminated with SpO2-style
    dropout spikes to ~0, must: (a) fail goodness-of-fit RAW, and
    (b) PASS with a significant trend after floor+MAD cleaning."""
    rng = np.random.default_rng(7)
    n = 4000
    spo2 = 94.0 + ci.simulate_approaching_fold(n, 1.0, 0.35, 0.03, 1.0, rng)
    drop = rng.choice(n, int(0.015 * n), replace=False)
    spo2[drop] = rng.uniform(0, 3, len(drop))

    raw = ci.ou_goodness_of_fit(ci.Series(np.arange(n, dtype=float), spo2, 1.0))
    assert not raw.adequate            # artifact must break the raw fit

    cleaned, n_removed = ci.clean_dropouts(spo2, floor=1.0, report=False)
    assert n_removed >= len(drop) // 2
    ser = ci.Series(np.arange(n, dtype=float), cleaned, 1.0)
    gof = ci.ou_goodness_of_fit(ser)
    tt = ci.surrogate_significance(ser, n_surrogates=200,
                                   rng=np.random.default_rng(0))
    assert gof.adequate                # cleaned fit must pass
    assert tt.tau_ar1 > 0.5 and tt.p_ar1 < 0.05   # signal recovered

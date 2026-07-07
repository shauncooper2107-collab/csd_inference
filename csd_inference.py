"""
csd_inference.py
================
An inference engine for Critical Slowing Down (CSD) and Early Warning
Signals (EWS) in real time series.

The tool inverts the usual simulation direction. Instead of choosing
parameters and watching a system evolve, it ingests an observed series
x(t) and *estimates* the dynamics, then asks the questions a complexity
scientist actually needs answered:

  1. What is the restoring rate k and noise level sigma of the process,
     with honest uncertainty?                         -> estimate_ou
  2. Do the standard CSD indicators (lag-1 autocorrelation, variance)
     trend upward as a transition is approached?      -> rolling_indicators
  3. Is that trend SIGNIFICANT, or could stationary noise produce it by
     chance?                                          -> surrogate_significance
  4. Does the linear (OU) framework even FIT this series, or is CSD the
     wrong lens for it?                                -> ou_goodness_of_fit
  5. Is the result robust to analyst choices (window length, detrending
     bandwidth)?                                      -> sensitivity_analysis
  6. Given the estimated dynamics, what is the survival probability of a
     corrective intervention applied over a finite window at a noise
     cost -- i.e. intervention as a first-passage problem?
                                                       -> intervention_forecast

Design commitments (why this is defensible):

  * The stochastic backbone is the Ornstein-Uhlenbeck process, the exact
    local linearisation of any smooth system near a stable fixed point:
        dX = -k (X - mu) dt + sigma dW
    Its exact discrete transition (Gaussian, conditional) gives a proper
    maximum-likelihood estimator for (k, mu, sigma), not a heuristic.

  * Significance is established against SURROGATES, not asymptotic
    hand-waving. A rising autocorrelation trend means nothing without a
    null that shares the series' spectrum but carries no trend.

  * The tool can say NO. Residual whiteness and normality tests report
    when the OU/CSD framework does not apply, so the engine does not
    manufacture a warning from every dataset.

  * Everything reproducible: a single seed governs all randomness.

This module is domain-general. It does not know or care whether the
series is a lake's turbidity, a patient's heart-rate variability, a
yeast population's density, or a market stress index. The organizational
front-end is a separate, explicitly-labelled interpretation layer and is
NOT part of this engine.

References the methods align with (for the reader to check against):
  Scheffer et al. (2009) Nature, "Early-warning signals for critical
  transitions"; Dakos et al. (2012) PLoS ONE, "Methods for detecting
  early warnings ..." and the associated `earlywarnings` R package.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy import stats
from scipy.ndimage import gaussian_filter1d


# ===========================================================================
#  data container
# ===========================================================================
@dataclass
class Series:
    """A uniformly-sampled scalar time series."""
    t: np.ndarray
    x: np.ndarray
    dt: float
    name: str = "series"

    @property
    def n(self) -> int:
        return len(self.x)


def load_series(path: str, time_col: Optional[str] = None,
<<<<<<< HEAD
                value_col: Optional[str] = None,
                dt_override: Optional[float] = None) -> Series:
    """Load a CSV. If column names are omitted, the first numeric column is
    taken as the value and (if a second exists) the first as time.

    dt_override sets the physical sampling interval explicitly and is the
    recommended path for real instrument data, where timestamps are often
    rounded, duplicated, or in awkward units. If it is omitted, dt is
    inferred from the time column and the loader REFUSES (rather than
    guessing) when that inference is ambiguous."""
=======
                value_col: Optional[str] = None) -> Series:
    """Load a CSV. If column names are omitted, the first numeric column is
    taken as the value and (if a second exists) the first as time. Uneven
    sampling is detected and reported; the analysis assumes uniform dt, so
    the engine warns rather than silently interpolating."""
>>>>>>> d69fdd40f49360b1245eb7faede3baa3635a57fc
    import pandas as pd
    df = pd.read_csv(path)
    numeric = df.select_dtypes(include=[np.number])
    if numeric.shape[1] == 0:
        raise ValueError("No numeric columns found in CSV.")
    if value_col is None:
        value_col = numeric.columns[-1]
    x = df[value_col].to_numpy(dtype=float)
    if time_col is not None:
        t = df[time_col].to_numpy(dtype=float)
    elif numeric.shape[1] >= 2:
        t = numeric.iloc[:, 0].to_numpy(dtype=float)
    else:
        t = np.arange(len(x), dtype=float)

    good = np.isfinite(x) & np.isfinite(t)
    x, t = x[good], t[good]
    order = np.argsort(t)
    t, x = t[order], x[order]

<<<<<<< HEAD
    # Determine the sampling interval. Silently defaulting dt=1.0 is dangerous:
    # every rate (k, tau), boundary, and forecast is expressed in units of dt,
    # so a wrong dt invalidates the whole report without any error. We therefore
    # REQUIRE an unambiguous dt: either supplied by the caller, or cleanly
    # inferable from the time column. Otherwise we refuse.
    if dt_override is not None:
        dt = float(dt_override)
        t = np.arange(len(x), dtype=float) * dt
    else:
        diffs = np.diff(t)
        med = float(np.median(diffs)) if diffs.size else 0.0
        if med <= 0:
            raise ValueError(
                "Cannot infer a sampling interval: the median timestamp "
                "difference is <= 0 (duplicate or rounded timestamps). "
                "Pass an explicit interval with --dt (e.g. --dt 1.0 for 1 Hz "
                "data), giving the real physical spacing between samples so "
                "that k, tau and the forecast are in meaningful units.")
        dt = med
        if np.std(diffs) > 0.05 * abs(dt):
            print(f"[warn] sampling looks uneven (dt spread "
                  f"{np.std(diffs):.3g} vs median {dt:.3g}). Analysis assumes "
                  f"uniform spacing at dt={dt:.4g}; consider resampling, or "
                  f"pass --dt to set the interval explicitly.")
    return Series(t=t, x=x, dt=dt, name=value_col)


# ===========================================================================
#  artifact / dropout cleaning
# ===========================================================================
def clean_dropouts(x: np.ndarray, floor: Optional[float] = None,
                   mad_threshold: float = 6.0,
                   report: bool = True) -> tuple:
    """Remove sensor-dropout and impulse artifacts before analysis.

    Real instrument channels (pulse oximetry, ECG-derived signals, etc.)
    produce non-physiological spikes when the probe detaches or the subject
    moves -- e.g. SpO2 plunging to 0 and recovering in a single sample. Left
    in, these dominate the residual variance and destroy the whiteness /
    normality tests, so the goodness-of-fit gate fires on the ARTIFACT rather
    than on the physiology.

    Two mechanisms:
      * floor: values at or below this are treated as 'no signal' and dropped
        (e.g. floor=1 for SpO2, since 0% is a sensor code, not a reading).
      * MAD spike test: points more than `mad_threshold` robust-SDs from a
        local median are flagged as impulses. MAD (median absolute deviation)
        is used instead of the SD because the SD is itself inflated by the
        spikes it is meant to catch.

    Flagged samples are linearly interpolated from their neighbours so the
    series stays uniformly sampled. Returns (cleaned_x, n_removed). This is
    reported transparently, never silent."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    flagged = np.zeros(n, dtype=bool)

    if floor is not None:
        flagged |= (x <= floor)

    # robust local spike detection on a median-filtered baseline
    from scipy.ndimage import median_filter
    win = max(5, min(51, (n // 100) | 1))       # odd window ~1% of series
    baseline = median_filter(x, size=win, mode="nearest")
    resid = x - baseline
    mad = np.median(np.abs(resid - np.median(resid)))
    if mad > 0:
        robust_sd = 1.4826 * mad                # MAD -> SD for Gaussian
        flagged |= np.abs(resid) > mad_threshold * robust_sd

    n_removed = int(flagged.sum())
    if n_removed == 0:
        return x.copy(), 0

    cleaned = x.copy()
    idx = np.arange(n)
    keep = ~flagged
    if keep.sum() < 2:
        raise ValueError("Nearly all samples flagged as artifacts; check "
                         "the floor / threshold settings.")
    cleaned[flagged] = np.interp(idx[flagged], idx[keep], x[keep])
    if report:
        print(f"[clean] removed {n_removed} artifact sample(s) "
              f"({100*n_removed/n:.2f}%) via "
              f"{'floor+' if floor is not None else ''}MAD spike test; "
              f"interpolated from neighbours.")
    return cleaned, n_removed
=======
    diffs = np.diff(t)
    dt = float(np.median(diffs))
    if diffs.size and np.std(diffs) > 0.05 * abs(dt):
        print(f"[warn] sampling looks uneven (dt spread "
              f"{np.std(diffs):.3g} vs median {dt:.3g}). Analysis assumes "
              f"uniform spacing; consider resampling.")
    return Series(t=t, x=x, dt=dt if dt != 0 else 1.0, name=value_col)
>>>>>>> d69fdd40f49360b1245eb7faede3baa3635a57fc


# ===========================================================================
#  detrending
# ===========================================================================
def detrend(x: np.ndarray, bandwidth: float) -> np.ndarray:
    """Gaussian-kernel detrending (the `earlywarnings` default idea).
    `bandwidth` is the kernel sigma in samples. Returns residuals. A slow
    trend must be removed or it inflates variance and biases autocorrelation."""
    if bandwidth <= 0:
        return x - np.mean(x)
    trend = gaussian_filter1d(x, sigma=bandwidth, mode="nearest")
    return x - trend


# ===========================================================================
#  OU parameter estimation (exact conditional MLE)
# ===========================================================================
@dataclass
class OUFit:
    k: float               # restoring rate (1/time). k -> 0 at transition.
    mu: float              # equilibrium level
    sigma: float           # noise intensity
    ar1: float             # e^{-k dt}, the lag-1 autocorrelation
    tau: float             # relaxation time 1/k
    stationary_sd: float   # sigma / sqrt(2k)
    k_ci: tuple = (np.nan, np.nan)
    sigma_ci: tuple = (np.nan, np.nan)
    n_used: int = 0


def _ou_point_estimate(x: np.ndarray, dt: float):
    """Exact-discretisation OU estimator. For dX=-k(X-mu)dt+sigma dW the
    exact transition is X_{i+1} = mu + (X_i-mu) B + eps, with
    B = e^{-k dt} and Var(eps) = sigma^2/(2k) (1 - B^2). Conditional least
    squares of X_{i+1} on X_i is the Gaussian MLE."""
    x0, x1 = x[:-1], x[1:]
    B, alpha, *_ = stats.linregress(x0, x1)
    # linregress returns slope, intercept, r, p, se -> slope=B, intercept=alpha
    B = float(B)
    alpha = float(alpha)
    if not (0.0 < B < 1.0):
        # No mean reversion detectable (B>=1) or anti-persistent (B<=0):
        # OU/CSD framework is questionable; flag with k at floor.
        B = min(max(B, 1e-6), 1 - 1e-6)
    k = -np.log(B) / dt
    mu = alpha / (1.0 - B)
    resid = x1 - (alpha + B * x0)
    v_eps = np.var(resid, ddof=2)
    sigma2 = 2.0 * k * v_eps / (1.0 - B * B)
    sigma = np.sqrt(max(sigma2, 0.0))
    return k, mu, sigma, B, resid


def estimate_ou(series: Series, n_boot: int = 500,
                rng: Optional[np.random.Generator] = None) -> OUFit:
    """Point estimate plus parametric-bootstrap 95% CIs for k and sigma.
    The bootstrap simulates OU replicates at the fitted parameters and
    refits, propagating estimation uncertainty honestly."""
    rng = rng or np.random.default_rng(0)
    x, dt = series.x, series.dt
    k, mu, sigma, B, _ = _ou_point_estimate(x, dt)

    ks, sigmas = [], []
    n = len(x)
    sd_eps = sigma * np.sqrt((1 - B * B) / (2 * k)) if k > 0 else 0.0
    for _ in range(n_boot):
        sim = np.empty(n)
        sim[0] = rng.normal(mu, sigma / np.sqrt(2 * k) if k > 0 else 1.0)
        noise = rng.normal(0.0, sd_eps, n - 1)
        for i in range(1, n):
            sim[i] = mu + (sim[i - 1] - mu) * B + noise[i - 1]
        kk, _, ss, *_ = _ou_point_estimate(sim, dt)
        ks.append(kk)
        sigmas.append(ss)
    k_ci = tuple(np.percentile(ks, [2.5, 97.5]))
    sigma_ci = tuple(np.percentile(sigmas, [2.5, 97.5]))

    return OUFit(k=k, mu=mu, sigma=sigma, ar1=B, tau=1.0 / k if k > 0 else np.inf,
                 stationary_sd=sigma / np.sqrt(2 * k) if k > 0 else np.inf,
                 k_ci=k_ci, sigma_ci=sigma_ci, n_used=n - 1)


# ===========================================================================
#  rolling indicators
# ===========================================================================
@dataclass
class RollingResult:
    centres: np.ndarray
    ar1: np.ndarray
    variance: np.ndarray
    skew: np.ndarray
    k_est: np.ndarray          # -ln(ar1)/dt, the CSD restoring rate
    window: int
    bandwidth: float


def _lag1_autocorr(w: np.ndarray) -> float:
    if len(w) < 3:
        return np.nan
    a, b = w[:-1], w[1:]
    sa, sb = a.std(), b.std()
    if sa < 1e-12 or sb < 1e-12:
        return np.nan
    return float(np.mean((a - a.mean()) * (b - b.mean())) / (sa * sb))


def rolling_indicators(series: Series, window_frac: float = 0.5,
                       bandwidth: float = None,
                       compute_skew: bool = True) -> RollingResult:
    """Compute the canonical CSD indicators in a sliding window over the
    detrended residuals. window_frac is the window length as a fraction of
    the series; bandwidth defaults to 10% of the series length.

    Rolling AR(1), variance and skew are computed via cumulative sums so the
    whole sweep is O(n), which keeps surrogate significance testing (which
    repeats this thousands of times) tractable."""
    x, dt, n = series.x, series.dt, series.n
    bw = bandwidth if bandwidth is not None else max(1.0, 0.1 * n)
    r = detrend(x, bw)
    W = max(10, int(round(window_frac * n)))
    W = min(W, n)
    m = W - 1
    n_win = n - W + 1

    def _cumwin(v, width):
        # sum of v over each sliding window of given width, length n-width+1
        c = np.concatenate(([0.0], np.cumsum(v)))
        return c[width:] + 0.0 - c[:len(c) - width]

    # variance (ddof=1) over full window W
    s1 = _cumwin(r, W)
    s2 = _cumwin(r * r, W)
    var = (s2 - s1 * s1 / W) / (W - 1)

    # lag-1 autocorr: within each window use a=r[start:start+W-1], b=r[start+1:start+W]
    prod = r[:-1] * r[1:]
    sq = r * r
    sum_a = _cumwin(r[:-1], m)                 # over first m of each window
    sum_b = _cumwin(r[1:], m)                  # over last m of each window
    sum_ab = _cumwin(prod, m)
    sum_a2 = _cumwin(sq[:-1], m)
    sum_b2 = _cumwin(sq[1:], m)
    # these arrays have length (n-1)-m+1 = n-W+1 = n_win  ✓
    cov = sum_ab - sum_a * sum_b / m
    va = sum_a2 - sum_a * sum_a / m
    vb = sum_b2 - sum_b * sum_b / m
    denom = np.sqrt(np.clip(va * vb, 1e-300, None))
    ar1 = np.where(denom > 1e-12, cov / denom, np.nan)

    # skew per window (optional; skipped in the surrogate hot path)
    centres = np.arange(n_win) + W // 2
    if compute_skew:
        skew = np.empty(n_win)
        for i in range(n_win):
            skew[i] = stats.skew(r[i:i + W])
    else:
        skew = np.full(n_win, np.nan)

    with np.errstate(divide="ignore", invalid="ignore"):
        kest = np.where((ar1 > 0) & (ar1 < 1), -np.log(ar1) / dt, np.nan)
    return RollingResult(centres=centres, ar1=ar1, variance=var, skew=skew,
                         k_est=kest, window=W, bandwidth=bw)


# ===========================================================================
#  trend + surrogate significance  (the rigor centrepiece)
# ===========================================================================
@dataclass
class TrendTest:
    tau_ar1: float
    tau_var: float
    p_ar1: float
    p_var: float
    n_surrogates: int
<<<<<<< HEAD
    null_method: str = "fourier"
=======
>>>>>>> d69fdd40f49360b1245eb7faede3baa3635a57fc
    null_ar1: np.ndarray = field(default=None, repr=False)


def _kendall(indicator: np.ndarray) -> float:
    idx = np.arange(len(indicator))
    good = np.isfinite(indicator)
    if good.sum() < 4:
        return np.nan
    tau, _ = stats.kendalltau(idx[good], indicator[good])
    return float(tau)


def fourier_surrogate(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Phase-randomised (Fourier) surrogate. Preserves the power spectrum --
    hence the linear autocorrelation structure -- while destroying any
<<<<<<< HEAD
    temporal TREND in that structure. This is the stricter, CSD-appropriate
    null: 'is the rise in autocorrelation real, or just what a STATIONARY
    correlated process produces by chance?'"""
=======
    temporal TREND in that structure. This is the correct null for 'is the
    rise in autocorrelation real, or just a stationary correlated process?'"""
>>>>>>> d69fdd40f49360b1245eb7faede3baa3635a57fc
    n = len(x)
    xm = x - x.mean()
    F = np.fft.rfft(xm)
    phases = rng.uniform(0, 2 * np.pi, size=F.shape)
    phases[0] = 0.0
    if n % 2 == 0:
        phases[-1] = 0.0
    F_s = np.abs(F) * np.exp(1j * phases)
    s = np.fft.irfft(F_s, n=n)
    return s + x.mean()


<<<<<<< HEAD
def permutation_surrogate(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Random permutation (shuffle) surrogate. Destroys ALL temporal
    structure, so the null is white noise. This is the conventional Dakos et
    al. (2012) test and the one most published warnings use; it is a WEAKER
    null than Fourier surrogates (any autocorrelation trend beats white
    noise), so it yields smaller p-values. Provided for direct comparability
    with the literature."""
    return rng.permutation(x)


_SURROGATE_FUNCS = {"fourier": fourier_surrogate,
                    "permutation": permutation_surrogate}


def surrogate_significance(series: Series, window_frac: float = 0.5,
                           bandwidth: float = None, n_surrogates: int = 1000,
                           null_method: str = "fourier",
                           rng: Optional[np.random.Generator] = None) -> TrendTest:
    """Kendall-tau trend of the rolling AR(1) and variance indicators, with
    one-sided p-values from surrogates. p = P(tau_null >= tau_obs).

    null_method selects the surrogate construction:
      'fourier'     -- stationary same-spectrum null (stricter; default).
                       Asks whether autocorrelation is *increasing* beyond a
                       stationary correlated process. Correct for CSD.
      'permutation' -- shuffle / white-noise null (conventional Dakos 2012).
                       Weaker; matches most published warnings. Use for
                       direct comparison with the literature.

    Report both when benchmarking: agreement is reassuring, and the gap is
    itself informative about how much of the 'signal' is merely correlation."""
    if null_method not in _SURROGATE_FUNCS:
        raise ValueError(f"null_method must be one of {list(_SURROGATE_FUNCS)}")
    surrogate = _SURROGATE_FUNCS[null_method]
=======
def surrogate_significance(series: Series, window_frac: float = 0.5,
                           bandwidth: float = None, n_surrogates: int = 1000,
                           rng: Optional[np.random.Generator] = None) -> TrendTest:
    """Kendall-tau trend of the rolling AR(1) and variance indicators, with
    one-sided p-values from Fourier surrogates. p = P(tau_null >= tau_obs).
    Small p => the upward CSD trend is unlikely under a stationary process
    with the same spectrum."""
>>>>>>> d69fdd40f49360b1245eb7faede3baa3635a57fc
    rng = rng or np.random.default_rng(0)
    obs = rolling_indicators(series, window_frac, bandwidth)
    tau_ar1 = _kendall(obs.ar1)
    tau_var = _kendall(obs.variance)

    null_ar1 = np.empty(n_surrogates)
    null_var = np.empty(n_surrogates)
    for s in range(n_surrogates):
<<<<<<< HEAD
        surr = Series(series.t, surrogate(series.x, rng),
=======
        surr = Series(series.t, fourier_surrogate(series.x, rng),
>>>>>>> d69fdd40f49360b1245eb7faede3baa3635a57fc
                      series.dt, series.name)
        r = rolling_indicators(surr, window_frac, bandwidth, compute_skew=False)
        null_ar1[s] = _kendall(r.ar1)
        null_var[s] = _kendall(r.variance)

    p_ar1 = (1 + np.sum(null_ar1 >= tau_ar1)) / (1 + n_surrogates)
    p_var = (1 + np.sum(null_var >= tau_var)) / (1 + n_surrogates)
    return TrendTest(tau_ar1=tau_ar1, tau_var=tau_var, p_ar1=p_ar1,
<<<<<<< HEAD
                     p_var=p_var, n_surrogates=n_surrogates,
                     null_method=null_method, null_ar1=null_ar1)
=======
                     p_var=p_var, n_surrogates=n_surrogates, null_ar1=null_ar1)
>>>>>>> d69fdd40f49360b1245eb7faede3baa3635a57fc


# ===========================================================================
#  goodness-of-fit  (the 'can say no' feature)
# ===========================================================================
@dataclass
class FitDiagnostics:
    ljung_box_p: float     # H0: residuals white. small p => OU misses structure
    normality_p: float     # H0: residuals Gaussian
    residual_ac1: float
    adequate: bool
    note: str


def ou_goodness_of_fit(series: Series, lags: int = 10) -> FitDiagnostics:
    """Fit OU, then test whether its residuals are white and Gaussian. If
    they are not, the linear/OU CSD framework is the wrong model for this
    series and any 'early warning' read from it is unreliable. A tool that
    cannot fail this test cannot be trusted when it passes."""
    x, dt = series.x, series.dt
    _, _, _, _, resid = _ou_point_estimate(x, dt)
    n = len(resid)

    # Ljung-Box Q on residual autocorrelation
    r = resid - resid.mean()
    denom = np.sum(r * r)
    acf = [np.sum(r[k:] * r[:-k]) / denom for k in range(1, lags + 1)]
    Q = n * (n + 2) * np.sum([(acf[k] ** 2) / (n - k - 1)
                              for k in range(lags)])
    lb_p = float(stats.chi2.sf(Q, df=lags))

    norm_p = float(stats.normaltest(resid).pvalue)
    adequate = (lb_p > 0.05) and (norm_p > 0.05)
    if adequate:
        note = "OU residuals consistent with white Gaussian noise: framework applies."
    elif lb_p <= 0.05:
        note = ("Residuals retain autocorrelation (Ljung-Box p<=0.05): the OU "
                "model misses structure (periodicity, higher-order dynamics). "
                "Treat CSD indicators with caution.")
    else:
        note = ("Residuals non-Gaussian (normaltest p<=0.05): heavy tails or "
                "jumps. Variance-based EWS may be unreliable; consider a "
                "jump-diffusion model.")
    return FitDiagnostics(ljung_box_p=lb_p, normality_p=norm_p,
                          residual_ac1=float(acf[0]), adequate=adequate,
                          note=note)


# ===========================================================================
#  sensitivity analysis
# ===========================================================================
@dataclass
class SensitivityResult:
    window_fracs: np.ndarray
    bandwidths: np.ndarray
    tau_grid: np.ndarray          # [len(window), len(bandwidth)]
    frac_positive: float
    robust: bool


def sensitivity_analysis(series: Series,
                         window_fracs=(0.25, 0.4, 0.5, 0.6, 0.75),
                         bandwidth_fracs=(0.05, 0.1, 0.2, 0.4)) -> SensitivityResult:
    """Recompute the AR(1) Kendall-tau across a grid of analyst choices.
    A trustworthy signal is stable across the grid; one that flips sign as
    the window nudges from 0.5 to 0.6 is an artefact. Reports the fraction
    of the grid with a positive trend."""
    n = series.n
    wf = np.array(window_fracs)
    bws = np.array([bf * n for bf in bandwidth_fracs])
    grid = np.full((len(wf), len(bws)), np.nan)
    for i, w in enumerate(wf):
        for j, bw in enumerate(bws):
            r = rolling_indicators(series, w, bw)
            grid[i, j] = _kendall(r.ar1)
    frac_pos = float(np.mean(grid > 0))
    return SensitivityResult(window_fracs=wf, bandwidths=bws, tau_grid=grid,
                             frac_positive=frac_pos, robust=(frac_pos > 0.9
                             or frac_pos < 0.1))


# ===========================================================================
#  intervention forecast  (the novel layer; bridges to the Kramers paper)
# ===========================================================================
def _mc_survival(k0, k1, sigma, d, H0, boundary, beta, T_ramp,
                 n_paths=400, dt=0.1, settle_frac=0.5, rng=None):
    """Monte-Carlo survival probability of a corrective intervention that
    ramps the restoring rate k0->k1 over window T while noise is transiently
    amplified to sigma*(1+beta*sin(pi t/T)). Truth for the time-inhomogeneous
    SDE; this is intervention framed as a first-passage problem."""
    rng = rng or np.random.default_rng(0)
    total = T_ramp * (1 + settle_frac)
    steps = max(int(total / dt), 2)
    sqdt = np.sqrt(dt)
    surv = 0
    for _ in range(n_paths):
        H = H0
        alive = True
        for s in range(steps):
            frac = min((s * dt) / T_ramp, 1.0)
            k = k0 + (k1 - k0) * frac
            spike = beta * np.sin(np.pi * frac) if frac < 1.0 else 0.0
            sig = sigma * (1 + spike)
            H += (-k * (H - 0.0) - d) * dt + sig * sqdt * rng.standard_normal()
            if H <= boundary:
                alive = False
                break
        if alive:
            surv += 1
    return surv / n_paths


@dataclass
class InterventionForecast:
    boundary: float
    p_survive: float
    p_ci: tuple
    k_from: float
    k_to: float
    beta: float
    T_ramp: float


def intervention_forecast(fit: OUFit, boundary_sds: float = 3.0,
                          k_multiplier: float = 2.0, beta: float = 1.0,
                          T_ramp: float = 20.0, drift: float = 0.0,
                          n_paths: int = 400,
                          rng: Optional[np.random.Generator] = None
                          ) -> InterventionForecast:
    """Using ESTIMATED dynamics, forecast the survival probability of a
    control that lifts k by k_multiplier over window T_ramp at noise cost
    beta, where 'collapse' is a boundary boundary_sds stationary-SDs below
    the equilibrium. Answers: can this system survive its own cure?"""
    rng = rng or np.random.default_rng(0)
    boundary = -boundary_sds * fit.stationary_sd
    p = _mc_survival(k0=fit.k, k1=fit.k * k_multiplier, sigma=fit.sigma,
                     d=drift, H0=0.0, boundary=boundary, beta=beta,
                     T_ramp=T_ramp, n_paths=n_paths, rng=rng)
    # Wilson 95% interval
    z, nn = 1.96, n_paths
    denom = 1 + z * z / nn
    centre = (p + z * z / (2 * nn)) / denom
    half = z * np.sqrt(p * (1 - p) / nn + z * z / (4 * nn * nn)) / denom
    return InterventionForecast(boundary=boundary, p_survive=p,
                                p_ci=(max(0, centre - half), min(1, centre + half)),
                                k_from=fit.k, k_to=fit.k * k_multiplier,
                                beta=beta, T_ramp=T_ramp)


# ===========================================================================
#  full report
# ===========================================================================
def run_analysis(series: Series, seed: int = 0, n_surrogates: int = 1000,
<<<<<<< HEAD
                 make_plot: Optional[str] = None,
                 both_nulls: bool = True) -> dict:
    rng = np.random.default_rng(seed)
    fit = estimate_ou(series, rng=rng)
    gof = ou_goodness_of_fit(series)
    trend = surrogate_significance(series, n_surrogates=n_surrogates,
                                   null_method="fourier", rng=rng)
    trend_perm = (surrogate_significance(series, n_surrogates=n_surrogates,
                                         null_method="permutation", rng=rng)
                  if both_nulls else None)
=======
                 make_plot: Optional[str] = None) -> dict:
    rng = np.random.default_rng(seed)
    fit = estimate_ou(series, rng=rng)
    gof = ou_goodness_of_fit(series)
    trend = surrogate_significance(series, n_surrogates=n_surrogates, rng=rng)
>>>>>>> d69fdd40f49360b1245eb7faede3baa3635a57fc
    sens = sensitivity_analysis(series)
    forecast = intervention_forecast(fit, rng=rng)

    print(f"\n{'='*70}\nCSD INFERENCE REPORT  —  {series.name}"
          f"   (n={series.n}, dt={series.dt:.4g}, seed={seed})\n{'='*70}")
    print(f"\n[1] ESTIMATED DYNAMICS (Ornstein-Uhlenbeck MLE)")
    print(f"    restoring rate k = {fit.k:.5g}  "
          f"95% CI [{fit.k_ci[0]:.4g}, {fit.k_ci[1]:.4g}]")
    print(f"    relaxation time  = {fit.tau:.4g}   (larger => nearer transition)")
    print(f"    noise sigma      = {fit.sigma:.4g}  "
          f"95% CI [{fit.sigma_ci[0]:.4g}, {fit.sigma_ci[1]:.4g}]")
    print(f"    equilibrium mu   = {fit.mu:.4g}   lag-1 AR = {fit.ar1:.4f}")

    print(f"\n[2] MODEL ADEQUACY (can this tool say 'no'?)")
    print(f"    Ljung-Box p = {gof.ljung_box_p:.3f}   "
          f"normality p = {gof.normality_p:.3f}")
    print(f"    {'PASS' if gof.adequate else 'CAUTION'}: {gof.note}")

<<<<<<< HEAD
    star = lambda p: "***" if p < 0.01 else "**" if p < 0.05 else "ns"
    print(f"\n[3] CSD TREND + SURROGATE SIGNIFICANCE "
          f"({trend.n_surrogates} surrogates)")
    print(f"    AR(1) Kendall tau = {trend.tau_ar1:+.3f}    "
          f"variance Kendall tau = {trend.tau_var:+.3f}")
    print(f"    p-values:                    AR(1)      variance")
    print(f"      Fourier null (strict):     {trend.p_ar1:.4f} {star(trend.p_ar1):>4}"
          f"   {trend.p_var:.4f} {star(trend.p_var):>4}")
    if trend_perm is not None:
        print(f"      permutation null (Dakos):  "
              f"{trend_perm.p_ar1:.4f} {star(trend_perm.p_ar1):>4}"
              f"   {trend_perm.p_var:.4f} {star(trend_perm.p_var):>4}")
        print(f"    (Fourier is the stricter, CSD-appropriate null; "
              f"permutation matches the\n     conventional literature test. "
              f"Divergence => 'signal' is partly mere autocorrelation.)")
=======
    print(f"\n[3] CSD TREND + SURROGATE SIGNIFICANCE "
          f"({trend.n_surrogates} Fourier surrogates)")
    star = lambda p: "***" if p < 0.01 else "**" if p < 0.05 else "ns"
    print(f"    AR(1) Kendall tau = {trend.tau_ar1:+.3f}   "
          f"p = {trend.p_ar1:.4f}  {star(trend.p_ar1)}")
    print(f"    variance   tau    = {trend.tau_var:+.3f}   "
          f"p = {trend.p_var:.4f}  {star(trend.p_var)}")
>>>>>>> d69fdd40f49360b1245eb7faede3baa3635a57fc

    print(f"\n[4] ROBUSTNESS (grid over window x bandwidth)")
    print(f"    fraction of grid with positive AR(1) trend = "
          f"{sens.frac_positive:.0%}   "
          f"{'ROBUST' if sens.robust else 'FRAGILE — result depends on choices'}")

    print(f"\n[5] INTERVENTION FORECAST (first-passage of the cure)")
    print(f"    boundary at {forecast.boundary:.3g} "
          f"(3 stationary-SD below eq.)")
    print(f"    lift k {forecast.k_from:.4g} -> {forecast.k_to:.4g} over "
          f"T={forecast.T_ramp:.0f}, noise cost beta={forecast.beta}")
    print(f"    P(survive the intervention) = {forecast.p_survive:.2f} "
          f"[{forecast.p_ci[0]:.2f}, {forecast.p_ci[1]:.2f}]")

    verdict = ("SIGNAL: significant, robust CSD trend"
               if (trend.p_ar1 < 0.05 and sens.frac_positive > 0.9
                   and gof.adequate)
               else "NO CLEAR SIGNAL / model caution — see notes above")
    print(f"\n{'='*70}\nVERDICT: {verdict}\n{'='*70}\n")
<<<<<<< HEAD
    print("(Verdict uses the strict Fourier null. The conventional "
          "permutation test may\n differ; see [3].)\n")

    results = dict(fit=fit, gof=gof, trend=trend, trend_perm=trend_perm,
                   sens=sens, forecast=forecast)
=======

    results = dict(fit=fit, gof=gof, trend=trend, sens=sens, forecast=forecast)
>>>>>>> d69fdd40f49360b1245eb7faede3baa3635a57fc
    if make_plot:
        _plot_report(series, results, make_plot)
    return results


def _plot_report(series: Series, res: dict, path: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fit, trend = res["fit"], res["trend"]
    obs = rolling_indicators(series)
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(f"CSD Inference — {series.name}", fontsize=13, weight="bold")

    ax[0, 0].plot(series.t, series.x, lw=0.8, color="#2c7fb8")
    trend_line = series.x - detrend(series.x, max(1.0, 0.1 * series.n))
    ax[0, 0].plot(series.t, trend_line, lw=1.5, color="#d95f0e",
                  label="Gaussian trend")
    ax[0, 0].set_title("Observed series + trend")
    ax[0, 0].legend(fontsize=8)

    c = obs.centres
    ax[0, 1].plot(c, obs.ar1, color="#e6550d")
    ax[0, 1].set_title(f"Rolling AR(1)   Kendall τ={trend.tau_ar1:+.2f}, "
                       f"p={trend.p_ar1:.3f}")
    ax[0, 1].set_ylabel("lag-1 autocorrelation")

    ax[1, 0].plot(c, obs.variance, color="#31a354")
    ax[1, 0].set_title(f"Rolling variance   Kendall τ={trend.tau_var:+.2f}, "
                       f"p={trend.p_var:.3f}")

    ax[1, 1].hist(trend.null_ar1, bins=40, color="#bdbdbd",
                  label="surrogate null")
    ax[1, 1].axvline(trend.tau_ar1, color="#e6550d", lw=2,
                     label=f"observed τ={trend.tau_ar1:+.2f}")
    ax[1, 1].set_title("AR(1) trend vs Fourier-surrogate null")
    ax[1, 1].legend(fontsize=8)

    for a in ax.flat:
        a.tick_params(labelsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"[plot] saved {path}")


# ===========================================================================
#  synthetic generators (for validation + demo)
# ===========================================================================
def simulate_ou(n, dt, k, sigma, mu=0.0, x0=0.0, rng=None):
    rng = rng or np.random.default_rng(0)
    B = np.exp(-k * dt)
    sd = sigma * np.sqrt((1 - B * B) / (2 * k))
    x = np.empty(n)
    x[0] = x0
    for i in range(1, n):
        x[i] = mu + (x[i - 1] - mu) * B + rng.normal(0, sd)
    return x


def simulate_approaching_fold(n, dt, k_start, k_end, sigma, rng=None):
    """Non-stationary OU whose restoring rate decays linearly k_start->k_end,
    the textbook 'approaching a transition' scenario. AR(1) should rise."""
    rng = rng or np.random.default_rng(0)
    ks = np.linspace(k_start, k_end, n)
    x = np.empty(n)
    x[0] = 0.0
    for i in range(1, n):
        B = np.exp(-ks[i] * dt)
        sd = sigma * np.sqrt((1 - B * B) / (2 * ks[i]))
        x[i] = x[i - 1] * B + rng.normal(0, sd)
    return x


# ===========================================================================
#  CLI
# ===========================================================================
def main():
    ap = argparse.ArgumentParser(
        description="CSD / early-warning-signal inference on a time series.")
    ap.add_argument("csv", nargs="?", help="CSV file; omit to run the demo.")
    ap.add_argument("--time-col", default=None)
    ap.add_argument("--value-col", default=None)
<<<<<<< HEAD
    ap.add_argument("--dt", type=float, default=None,
                    help="physical sampling interval (e.g. 1.0 for 1 Hz). "
                         "Strongly recommended for real instrument data; "
                         "required when timestamps are rounded/duplicated.")
    ap.add_argument("--clean", action="store_true",
                    help="remove sensor-dropout / impulse artifacts before "
                         "analysis (MAD spike test).")
    ap.add_argument("--floor", type=float, default=None,
                    help="values <= floor are treated as 'no signal' and "
                         "cleaned (e.g. --floor 1 for SpO2). Implies --clean.")
=======
>>>>>>> d69fdd40f49360b1245eb7faede3baa3635a57fc
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--surrogates", type=int, default=1000)
    ap.add_argument("--plot", default=None, help="path to save a PNG report")
    args = ap.parse_args()

    if args.csv:
<<<<<<< HEAD
        series = load_series(args.csv, args.time_col, args.value_col,
                             dt_override=args.dt)
        if args.clean or args.floor is not None:
            cleaned, _ = clean_dropouts(series.x, floor=args.floor)
            series = Series(series.t, cleaned, series.dt, series.name)
=======
        series = load_series(args.csv, args.time_col, args.value_col)
>>>>>>> d69fdd40f49360b1245eb7faede3baa3635a57fc
        run_analysis(series, seed=args.seed, n_surrogates=args.surrogates,
                     make_plot=args.plot)
    else:
        print("No CSV given — running built-in validation demo.")
        _demo()


def _demo():
    """Three-part validation the reader can run with no data:
      (1) a stationary control the tool should NOT flag,
      (2) an approaching-fold series it SHOULD flag, and
      (3) a non-OU (periodic) series whose goodness-of-fit it should REJECT.
    A tool that passes all three earns trust on real data."""
    print("Built-in validation demo — three synthetic cases.\n")

    print("### CASE 1: stationary control (expect: NO signal) ###")
    x = simulate_ou(800, 1.0, k=0.2, sigma=1.0,
                    rng=np.random.default_rng(2002))
    run_analysis(Series(np.arange(800.0), x, 1.0, "CONTROL-stationary"),
                 seed=2, n_surrogates=500)

    print("\n### CASE 2: approaching a fold (expect: significant, robust) ###")
    fold = simulate_approaching_fold(800, 1.0, 0.4, 0.02, 1.0,
                                     np.random.default_rng(1))
    run_analysis(Series(np.arange(800.0), fold, 1.0, "APPROACHING-fold"),
                 seed=1, n_surrogates=500, make_plot="demo_report.png")

    print("\n### CASE 3: non-OU periodic series (expect: model REJECTED) ###")
    t = np.arange(800.0)
    osc = simulate_ou(800, 1.0, k=0.3, sigma=0.5,
                      rng=np.random.default_rng(4)) + 3 * np.sin(2 * np.pi * t / 12)
    run_analysis(Series(t, osc, 1.0, "NON-OU-periodic"),
                 seed=4, n_surrogates=500)


if __name__ == "__main__":
    main()

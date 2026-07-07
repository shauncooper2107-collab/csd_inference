# CSD Inference Engine

**Detecting critical slowing down and early-warning signals in real time series — with the statistical rigor needed to believe the result, and the honesty to say when there isn't one.**

This is a domain-general inference tool. Point it at a single observable measured through time — a lake's turbidity, a patient's inter-beat intervals, a yeast population's density, a market stress index — and it estimates the underlying dynamics and tests, rigorously, whether the system is losing resilience.

It does **not** simulate a metaphor. It **infers** from your data.

---

## Why this exists

Most early-warning-signal (EWS) work runs one direction: pick parameters, simulate, observe. That is useful for teaching and exploration but it cannot tell a scientist anything about *their* system. This tool inverts that: it takes an observed series and estimates the restoring rate `k`, the noise level `sigma`, and the equilibrium, then asks the questions that actually matter for a real dataset.

The design is built around one principle: **a tool that cannot fail is worthless when it passes.** Every result comes with a null model, an uncertainty, a robustness check, and a goodness-of-fit test that can reject the framework outright.

---

## What it computes

| Step | Question answered | Method |
|---|---|---|
| **Estimate dynamics** | What are `k`, `sigma`, equilibrium? | Exact-discretisation Ornstein–Uhlenbeck maximum likelihood, with parametric-bootstrap 95% CIs |
| **Model adequacy** | Does the OU/CSD framework even *fit* this series? | Ljung–Box test (residual whiteness) + normality test. **Can return CAUTION / reject.** |
| **CSD trend** | Do autocorrelation and variance rise? | Rolling lag-1 autocorrelation and variance; Kendall's τ trend |
| **Significance** | Is that rise beyond chance? | One-sided p-value from **Fourier phase-randomised surrogates** that share the spectrum but carry no trend |
| **Robustness** | Does the result survive analyst choices? | Kendall's τ recomputed over a grid of window lengths × detrending bandwidths |
| **Intervention forecast** | Can the system survive its own cure? | Monte-Carlo first-passage survival of a control that lifts `k` over a finite window at a noise cost |

The last row is the novel contribution: interventions are modelled as a **first-passage problem during a parameter ramp**. Standard EWS toolkits tell you a transition is coming; none model the survival probability of the correction itself.

---

## Install

```
pip install numpy scipy matplotlib pandas
```

That's the whole dependency list — the standard scientific stack.

---

## Use it on your data

Your CSV needs one column of measurements, optionally a time column:

```
time,value
0.0,12.4
1.0,12.9
2.0,11.8
...
```

Then:

```
python csd_inference.py mydata.csv --plot report.png
```

or specify columns explicitly:

```
python csd_inference.py mydata.csv --value-col turbidity --time-col day --plot report.png
```

You get a printed report (dynamics, adequacy, trend + significance, robustness, intervention forecast, and a plain verdict) plus a four-panel PNG: the series with its trend, rolling AR(1), rolling variance, and your observed τ against the surrogate null.

**Reproducibility:** every run is governed by a single `--seed`. The same seed reproduces every number, including the Monte-Carlo forecast and the surrogate p-values.

### From Python

```python
import csd_inference as ci

series  = ci.load_series("mydata.csv")
fit     = ci.estimate_ou(series)              # k, sigma, CIs
gof     = ci.ou_goodness_of_fit(series)       # does the model fit?
trend   = ci.surrogate_significance(series)   # tau + p-value
sens    = ci.sensitivity_analysis(series)     # robustness grid
cast    = ci.intervention_forecast(fit)       # survival of a correction

results = ci.run_analysis(series, make_plot="report.png")  # all of it
```

---

## Real instrument data: sampling interval and artifacts

Real monitoring channels (pulse oximetry, ECG-derived signals, sensor feeds) bring two problems that will silently wreck an EWS analysis if ignored. The tool now guards against both.

### Set the sampling interval explicitly

Every rate the engine reports — `k`, the relaxation time `τ`, the intervention window `T`, the collapse boundary — is expressed in units of the sampling interval `dt`. If `dt` is wrong, all of them are wrong, with no error to warn you.

Instrument timestamps are often rounded or duplicated, which makes `dt` impossible to infer reliably (the median gap can collapse to zero). **The loader now refuses to guess** in that case and asks you to state the real interval:

```
python csd_inference.py signal.csv --dt 1.0     # 1 Hz data
python csd_inference.py signal.csv --dt 0.004   # 250 Hz data
```

Pass the true physical spacing between samples. When in doubt, `--dt` is always safer than letting the timestamps decide.

### Clean sensor dropouts before analysis

A probe detaching or a patient moving produces non-physiological spikes — SpO2 plunging to 0% and recovering in a single sample, for instance. Left in, these dominate the residual variance and trip the goodness-of-fit gate on the *artifact* rather than the physiology, burying any real signal.

```
python csd_inference.py spo2.csv --dt 1.0 --floor 1     # 0% = sensor code, not a reading
python csd_inference.py ecg.csv  --dt 0.004 --clean      # MAD spike test only
```

- `--floor V` treats values at or below `V` as "no signal" and removes them.
- `--clean` runs a robust MAD (median-absolute-deviation) spike test that flags impulses without being fooled by the spikes it is trying to catch.

Flagged samples are linearly interpolated from their neighbours, and the number removed is always reported — never silent.

**Honest caveat:** interpolating across gaps can itself smooth the series and inflate autocorrelation. Treat a signal that appears only after cleaning as *worth investigating*, not proven. Confirm it does not depend on the interpolation — for example by comparing against simply excising the artifact segments.

---

Run the engine with no arguments to execute a built-in three-case validation:

```
python csd_inference.py
```

1. **Stationary control** — a constant-`k` process. The tool should recover `k` and return **no significant trend**. (It does: τ ≈ 0.23, p ≈ 0.39, ns.)
2. **Approaching a fold** — `k` decays toward zero. The tool should return a **significant, robust** signal. (It does: τ ≈ 0.97, p ≈ 0.002, robust across the whole grid.)
3. **Non-OU periodic series** — dynamics the OU model cannot capture. The tool should **reject its own framework**. (It does: Ljung–Box p < 0.0001 → CAUTION.)

Calibration was checked directly: on independent stationary series the false-positive rate is at or below the nominal α (≈ 2.5% at α = 0.05, p-values approximately uniform). The tool does not manufacture warnings.

---

## Two significance nulls

Significance is established against surrogates, and the tool reports two:

- **Fourier null (default, strict):** phase-randomised surrogates that preserve the series' spectrum (its autocorrelation) but destroy any *trend* in it. This asks the CSD-appropriate question — is autocorrelation *increasing* beyond what a stationary correlated process produces? The verdict uses this null.
- **Permutation null (conventional):** shuffled surrogates (white noise), the standard Dakos et al. (2012) test that most published warnings use. Weaker, so it gives smaller p-values.

Both are printed side by side. Agreement is reassuring; a gap tells you how much of the apparent "signal" is merely the series' own autocorrelation rather than a genuine rising trend. Select one in code with `null_method="fourier"` or `"permutation"`.

## Benchmark against the literature

The `benchmarks/` folder holds real reference series from the critical-transitions literature (converted from the `earlywarnings` R package):

```
python csd_inference.py benchmarks/foldbif.csv --dt 1.0 --surrogates 500
```

On `foldbif` (the canonical fold bifurcation) the engine recovers AR(1) Kendall τ ≈ +0.72, matching the published value — validating the indicator computation against a peer-reviewed fixture. Significance is marginal on this short (n=200) series under both nulls, which is honest rather than tuned. See `benchmarks/README.md`.

---

## What it deliberately does *not* do

- It does not interpret your series for you. `k` is a restoring rate, not a diagnosis.
- It does not claim power it lacks. EWS detection is genuinely hard; on short or noisy series the honest answer is often "no clear signal," and the tool will say so.
- It contains **no organizational, financial, or ecological interpretation layer.** Those belong in a separate, explicitly-labelled application built *on top of* this engine, never inside it.

---

## Method notes for reviewers

- **OU estimator.** For `dX = -k(X-μ)dt + σ dW`, the exact discrete transition is `X_{i+1} = μ + (X_i-μ)e^{-kΔt} + ε`, `Var(ε) = σ²/(2k)(1-e^{-2kΔt})`. Conditional least squares of `X_{i+1}` on `X_i` is the Gaussian MLE; `k = -ln(slope)/Δt`. CIs by parametric bootstrap.
- **Surrogates.** Fourier phase-randomisation preserves the periodogram (hence the linear autocorrelation) while destroying any temporal *trend* in autocorrelation — the correct null for "is the rise in AR(1) real, or a stationary correlated process?"
- **Detrending.** Gaussian-kernel, bandwidth reported and varied in the sensitivity grid.
- **Intervention forecast.** Time-inhomogeneous SDE integrated by Euler–Maruyama; survival = fraction of paths not crossing an absorbing boundary set a chosen number of stationary SDs below equilibrium; Wilson interval on the proportion. The quasi-static mean-first-passage-time approximation (exact OU escape integral) is available as a theory comparison.

The alignment target is the standard literature: Scheffer et al. (2009) for the indicator set and Dakos et al. (2012) / the `earlywarnings` package for detrending, rolling indicators, Kendall-τ trend, and sensitivity analysis. This engine's addition to that toolkit is the intervention first-passage layer and the built-in goodness-of-fit gate.

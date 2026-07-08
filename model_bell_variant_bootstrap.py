#!/usr/bin/env python3
"""
model_bell_variant_bootstrap.py
================================

Bootstrap / leave-one-out cross-validation robustness check for the bell-
shaped (Gaussian-in-S) light-field variant's fitted peak location, S_peak.

Builds ON TOP OF model_bell_variant.py (imported here, NOT modified or
duplicated) -- see that script for the full physics writeup, citations, and
the single-fit result (S_peak=1.4, from a 5-start multi-start check) that
this script stress-tests with proper resampling statistics.

WHY THIS SCRIPT EXISTS: model_bell_variant.py's own calibrate-glycine mode
already flagged that its 4-parameter fit (kappa_max, S_peak, width, log10(A))
against only 3 real data points has NEGATIVE residual degrees of freedom (3-4
= -1), and its multi-start check found S_peak was NOT consistently identified
across 5 different initial guesses. This script goes further with an actual
bootstrap (>=500 resamples) and leave-one-out cross-validation, to turn
"different starts disagree" into a real sampling distribution for S_peak.

HEADLINE CAVEAT -- read this before the numbers below, not after:
The real calibration dataset has only N=3 points (S=1.4, 1.5, 1.6). A
standard case-resampling bootstrap draws N=3 points WITH REPLACEMENT from
those same 3 points, which means a large fraction of resamples will contain
only 1 or 2 UNIQUE supersaturation values, not 3. This script reports that
unique-point-count distribution explicitly, because it is the single biggest
reason to expect (and correctly interpret) whatever spread we see in S_peak
below. Bootstrapping N=3 data is already at (arguably below) the edge of
where nonparametric bootstrap theory is normally considered meaningful --
this script is being run anyway, precisely to quantify how bad the situation
actually is, not to pretend the bootstrap solves the small-N problem.
"""

import sys
import importlib.util
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent


def _load_module(filename, module_name):
    """Load a sibling script as an importable module without running its
    `if __name__ == '__main__'` block."""
    spec = importlib.util.spec_from_file_location(module_name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    # Must register in sys.modules BEFORE exec_module: @dataclass looks up its
    # own module via sys.modules to resolve type hints.
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


bell = _load_module("model_bell_variant.py", "nplin_bell_base")


# =============================================================================
# Shared setup: identical material/light/geometry placeholders to the ones
# model_bell_variant.py's own run_calibrate_glycine() uses, so this script is
# testing the SAME fit, not a different one.
# =============================================================================
def build_template_objects():
    ds = bell.GLYCINE_JAVID2016
    laser = ds["laser"]
    M_glycine = 75.07e-3   # kg/mol -- real, textbook value
    rho_glycine = 1.6e3    # kg/m^3 -- real, approximate literature value for glycine crystal density
    v_m = M_glycine / (rho_glycine * bell.sc.N_A)
    material = bell.MaterialParams(
        name="glycine, bell variant (bootstrap analysis)",
        gamma=1.3e-2, v_m=v_m, T=ds["T_K"], A_prefactor=1e20,
        n_solvent=1.33, n_solid=1.55,
    )
    V_assumed = bell.estimate_illuminated_volume_m3(ds)
    light_template = bell.LightFieldParams(
        wavelength_nm=laser["wavelength_nm"],
        pulse_energy_J=laser["power_density_W_m2"] * np.pi * (laser["spot_diameter_m"] / 2) ** 2 * laser["pulse_duration_s"],
        spot_diameter_m=laser["spot_diameter_m"], pulse_duration_s=laser["pulse_duration_s"],
        n_pulses=int(laser["exposure_time_s"] * laser["rep_rate_Hz"]),
        coupling_efficiency=1.0, pulse_shape_factor=1.0, S_peak=1.5, bell_width=0.1,
    )
    geom = bell.SampleGeometry(volume_m3=V_assumed, exposure_time_s=laser["exposure_time_s"])
    I_peak = laser["power_density_W_m2"]
    return ds, material, light_template, geom, I_peak


def get_real_data(ds):
    S_data, A_data, A_sigma = [], [], []
    for S, row in sorted(ds["table2_biexponential_fit"].items()):
        S_data.append(S)
        A_data.append(row["A"])
        A_sigma.append(row["A_ci95"] / 1.96)  # 95% CI half-width -> 1-sigma
    return np.array(S_data), np.array(A_data), np.array(A_sigma)


# =============================================================================
# Core single-fit routine (silent), reused for every bootstrap resample and
# every leave-one-out refit -- SAME fitting procedure and SAME plausibility
# bound as model_bell_variant.py's run_calibrate_glycine(), just without its
# per-run console narration (which would be unreadable across 500+ calls).
# =============================================================================
BOUNDS = [(0.0, 1.0e7), (1.0, 2.0), (0.02, 1.0), (0.0, 40.0)]
INITIAL_S_PEAKS = [1.2, 1.5, 1.8]  # fewer starts than the original script's 5, for bootstrap speed

# Coarser than model_bell_variant.py's 60-point constraint grid, purely for
# computational speed across hundreds of refits. The qualitative conclusions
# below do not depend sensitively on this resolution.
N_CONSTRAINT_GRID = 20


def model_bell_predict(S_arr, kappa_max, S_peak, width, log_A, material, geom, light_template,
                        I_peak, n_solvent, n_solid):
    A_pref = 10 ** log_A
    preds = []
    for S in S_arr:
        u_field = bell.u_field_bell(S, kappa_max, S_peak, width, I_peak, n_solvent, n_solid)
        J = bell.nucleation_rate_field(S, material.T, material.gamma, material.v_m, A_pref, u_field)
        preds.append(bell.probability_multi_pulse(J, geom.volume_m3, light_template.pulse_duration_s,
                                                    light_template.n_pulses))
    return np.array(preds)


def fit_bell_once(S_data, A_data, A_sigma, material, geom, light_template, I_peak, n_solvent,
                   n_solid, S_check_grid):
    """Multi-start constrained fit of (kappa_max, S_peak, width, log_A) to the given
    (possibly resampled/reduced) dataset, under the SAME plausibility-bound constraint
    used in model_bell_variant.py. Returns the best (lowest chi^2) converged result among
    the tried starting points, or None if every start failed to converge."""

    def objective_chi2(params):
        kappa_max, S_peak, width, log_A = params
        preds = model_bell_predict(S_data, kappa_max, S_peak, width, log_A, material, geom,
                                    light_template, I_peak, n_solvent, n_solid)
        return np.sum(((preds - A_data) / A_sigma) ** 2)

    def ratio_constraint(params):
        kappa_max, S_peak, width, log_A = params
        return bell.MAX_PLAUSIBLE_NPLIN_RATIO - bell.max_ratio_over_grid_bell(
            kappa_max, S_peak, width, S_check_grid, material, I_peak, n_solvent, n_solid)

    constraints = [{"type": "ineq", "fun": ratio_constraint}]
    results = []
    for S_peak0 in INITIAL_S_PEAKS:
        x0 = [10.0, S_peak0, 0.15, 20.0]
        res = minimize(objective_chi2, x0, method="SLSQP", bounds=BOUNDS,
                        constraints=constraints, options={"maxiter": 300, "ftol": 1e-10})
        if res.success:
            results.append(res)
    if not results:
        return None
    best = min(results, key=lambda r: r.fun)
    kappa_max, S_peak, width, log_A = best.x
    return {"kappa_max": kappa_max, "S_peak": S_peak, "width": width, "logA": log_A, "chi2": best.fun}


def kappa_at_bound(kappa_max):
    lo, hi = BOUNDS[0]
    return kappa_max < 1e-3 or kappa_max > hi * 0.999


def s_peak_at_bound(S_peak):
    lo, hi = BOUNDS[1]
    return (S_peak - lo) < 1e-3 or (hi - S_peak) < 1e-3


def width_at_bound(width):
    lo, hi = BOUNDS[2]
    return (width - lo) < 1e-3 or (hi - width) < 1e-3


def any_at_bound(result):
    return (kappa_at_bound(result["kappa_max"]) or s_peak_at_bound(result["S_peak"])
            or width_at_bound(result["width"]))


# =============================================================================
# Main analysis
# =============================================================================
def main():
    print("=" * 78)
    print("BOOTSTRAP / CROSS-VALIDATION ROBUSTNESS CHECK -- bell-shaped variant's S_peak")
    print("Builds on model_bell_variant.py (imported, not modified). See that script for")
    print("the full physics writeup and the single-fit result being stress-tested here.")
    print("=" * 78)

    ds, material, light_template, geom, I_peak = build_template_objects()
    S_data_full, A_data_full, A_sigma_full = get_real_data(ds)
    n_solvent, n_solid = material.n_solvent, material.n_solid
    S_check_grid = np.linspace(*bell.RATIO_CONSTRAINT_S_RANGE, N_CONSTRAINT_GRID)

    N = len(S_data_full)
    print(f"\n--- Dataset size: N = {N} real calibration points ---")
    print(f"S values: {S_data_full.tolist()}")
    print("THIS NUMBER MATTERS A LOT: a case-resampling bootstrap of N=3 draws 3 points WITH")
    print("REPLACEMENT from these same 3 -- many resamples will contain only 1 or 2 UNIQUE S")
    print("values, not 3. That fact, not any subtlety of the optimizer, is expected to be the")
    print("single biggest driver of whatever spread we see in S_peak below.")

    # ------------------------------------------------------------------
    # Bootstrap
    # ------------------------------------------------------------------
    N_BOOT = 500
    rng = np.random.default_rng(0)
    boot_results = []
    n_unique_counts = []
    n_failed = 0

    print(f"\n--- Running {N_BOOT} bootstrap resamples (may take a minute) ---")
    for b in range(N_BOOT):
        idx = rng.integers(0, N, size=N)
        n_unique_counts.append(len(set(idx.tolist())))
        S_b, A_b, sig_b = S_data_full[idx], A_data_full[idx], A_sigma_full[idx]
        result = fit_bell_once(S_b, A_b, sig_b, material, geom, light_template, I_peak,
                                n_solvent, n_solid, S_check_grid)
        if result is None:
            n_failed += 1
            continue
        result["at_bound"] = any_at_bound(result)
        result["n_unique"] = n_unique_counts[-1]
        boot_results.append(result)
        if (b + 1) % 100 == 0:
            print(f"  ...{b + 1}/{N_BOOT} resamples done")

    n_at_bound = sum(1 for r in boot_results if r["at_bound"])
    print(f"\nBootstrap complete: {len(boot_results)}/{N_BOOT} fits converged "
          f"({n_failed} failed to converge -- {100 * n_failed / N_BOOT:.1f}%)")
    print(f"Of the converged fits, {n_at_bound}/{max(len(boot_results), 1)} "
          f"({100 * n_at_bound / max(len(boot_results), 1):.1f}%) landed AT a parameter bound "
          f"(kappa_max, S_peak, or width pinned to its box limit -- not an interior optimum)")

    unique_counts_arr = np.array(n_unique_counts)
    print("\nUnique-S-value count across the 500 resamples (see headline caveat above):")
    for k in range(1, N + 1):
        frac = np.mean(unique_counts_arr == k)
        print(f"  resamples with exactly {k} unique S value(s): {100 * frac:.1f}%")

    if not boot_results:
        print("\nNO bootstrap fits converged at all. Cannot proceed to distribution/plot/verdict --")
        print("this alone is the strongest possible version of 'underdetermined by the data'.")
        return

    S_peaks = np.array([r["S_peak"] for r in boot_results])
    widths = np.array([r["width"] for r in boot_results])
    kappas = np.array([r["kappa_max"] for r in boot_results])

    def summarize(name, arr):
        mean, std = np.mean(arr), np.std(arr)
        p5, p95 = np.percentile(arr, [5, 95])
        print(f"{name}: mean={mean:.4g}, std={std:.4g}, 90% interval=[{p5:.4g}, {p95:.4g}]")
        return mean, std, p5, p95

    print("\n--- Bootstrap distribution summary (converged fits only) ---")
    S_peak_mean, S_peak_std, S_peak_p5, S_peak_p95 = summarize("S_peak", S_peaks)
    width_mean, width_std, width_p5, width_p95 = summarize("width", widths)
    kappa_mean, kappa_std, kappa_p5, kappa_p95 = summarize("kappa_max", kappas)

    # ------------------------------------------------------------------
    # Leave-one-out cross-validation
    # ------------------------------------------------------------------
    print(f"\n--- Leave-one-out cross-validation (N={N}, so {N} refits with {N - 1} points each) ---")
    print(f"NOTE: {N - 1} real data points against 4 free parameters is EVEN MORE underdetermined")
    print(f"than the full N={N} case (residual DOF = {N - 1} - 4 = {N - 1 - 4}) -- expect this to be")
    print("at least as degenerate as the bootstrap results, likely more so.")
    loo_results = []
    for i in range(N):
        mask = np.arange(N) != i
        result = fit_bell_once(S_data_full[mask], A_data_full[mask], A_sigma_full[mask],
                                material, geom, light_template, I_peak, n_solvent, n_solid,
                                S_check_grid)
        left_out = S_data_full[i]
        if result is None:
            print(f"  left out S={left_out:.2f}  ->  FAILED to converge")
            loo_results.append(None)
            continue
        flag = "  <- AT a parameter bound" if any_at_bound(result) else ""
        print(f"  left out S={left_out:.2f}  ->  S_peak={result['S_peak']:.4g}, "
              f"width={result['width']:.4g}, kappa_max={result['kappa_max']:.4g}{flag}")
        loo_results.append(result)

    # ------------------------------------------------------------------
    # Plot: histogram of bootstrapped S_peak, width, kappa_max
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    axes[0].hist(S_peaks, bins=30, color="tab:blue", alpha=0.75, edgecolor="k")
    axes[0].axvline(1.4, color="tab:purple", ls="-", lw=2.5, label="original single fit, S_peak=1.4")
    axes[0].axvspan(1.45, 1.55, color="tab:red", alpha=0.2, label="literature window S=1.45-1.55")
    axes[0].set_xlabel("Bootstrapped S_peak")
    axes[0].set_ylabel("Count")
    axes[0].set_title(f"S_peak (std={S_peak_std:.3g})", fontsize=10)
    axes[0].legend(fontsize=7, loc="upper right")

    axes[1].hist(widths, bins=30, color="tab:green", alpha=0.75, edgecolor="k")
    axes[1].axvline(0.02, color="tab:purple", ls="-", lw=2.5, label="original single fit, width=0.02")
    axes[1].set_xlabel("Bootstrapped width")
    axes[1].set_title(f"width (std={width_std:.3g})", fontsize=10)
    axes[1].legend(fontsize=7, loc="upper right")

    # Some bootstrap fits converge to kappa_max=0 exactly (i.e. "no light effect at all" is
    # the best fit for that resample) -- log10(0) = -inf, which breaks a log-scale histogram.
    # This is itself a meaningful finding, not just a plotting nuisance, so report it
    # explicitly rather than silently dropping it or letting it crash.
    n_kappa_zero = int(np.sum(kappas <= 1e-3))
    kappas_nonzero = kappas[kappas > 1e-3]
    print(f"\n{n_kappa_zero}/{len(boot_results)} bootstrap fits converged to kappa_max~0 "
          f"(no light-field effect at all was the best fit for that resample).")
    if len(kappas_nonzero) > 0:
        axes[2].hist(np.log10(kappas_nonzero), bins=30, color="tab:orange", alpha=0.75, edgecolor="k",
                     label=f"kappa_max>0 (n={len(kappas_nonzero)})")
    axes[2].axvline(np.log10(15.03), color="tab:purple", ls="-", lw=2.5,
                    label="original single fit, kappa_max=15.03")
    axes[2].set_xlabel("log10(bootstrapped kappa_max)")
    axes[2].set_title(f"kappa_max (log10 scale); {n_kappa_zero} fits had kappa_max=0", fontsize=9)
    axes[2].legend(fontsize=7, loc="upper right")

    fig.suptitle(f"Bootstrap distributions, N_real={N}, {N_BOOT} resamples "
                 f"({len(boot_results)} converged, {n_failed} failed, {n_at_bound} at a bound)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig("bell_bootstrap_histograms.png", dpi=150)
    print("\nSaved plot to bell_bootstrap_histograms.png")

    # ------------------------------------------------------------------
    # Explicit, unhedged verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)
    convergence_rate = len(boot_results) / N_BOOT
    at_bound_rate = n_at_bound / max(len(boot_results), 1)

    # Thresholds, stated explicitly rather than left implicit:
    #  - "most fits fail/hit bounds": >50% failure OR >50% of converged fits at a bound.
    #  - "tight clustering": std well below the ~0.05-0.15 gap between S_peak=1.4 and the
    #    literature window -- we use 0.05 as the cutoff (roughly the gap to the window's
    #    near edge).
    #  - "wide scatter": std comparable to a large fraction of the full [1.0, 2.0] search
    #    range -- we use 0.2 (a fifth of that range) as the cutoff.
    if convergence_rate < 0.5 or at_bound_rate > 0.5:
        print(f"MOST bootstrap fits are unreliable: {100 * (1 - convergence_rate):.0f}% failed to")
        print(f"converge, and {100 * at_bound_rate:.0f}% of those that DID converge pinned to a")
        print("parameter bound rather than an interior optimum.")
        print("\nVERDICT: The bell-shaped model is UNDERDETERMINED by the available data. A")
        print("confident peak location CANNOT currently be claimed at all -- not near S=1.4, not")
        print("near the literature window, not anywhere. This is a direct, quantified consequence")
        print("of having only 3 real calibration points for a 4-parameter model.")
    elif S_peak_std >= 0.2:
        print(f"S_peak is WIDELY SCATTERED across bootstrap resamples: std={S_peak_std:.3g}, "
              f"comparable to a large fraction of the entire S_peak search range [1.0, 2.0].")
        print("\nVERDICT: The original single fit's peak location (S_peak=1.4) is NOT reliable. It")
        print("was likely a product of limited data and an underdetermined optimization landing")
        print("wherever it started, not a real physical signal. Do not treat S_peak=1.4 -- or its")
        print("proximity or non-proximity to the literature window -- as a finding.")
    elif S_peak_std < 0.05:
        print(f"S_peak clusters relatively tightly across bootstrap resamples: std={S_peak_std:.3g},")
        print("clearly smaller than the ~0.05-0.15 gap between the original fit (S_peak=1.4) and")
        print("the literature window (S=1.45-1.55).")
        print("\nVERDICT: The original fit's peak location looks ROBUST to resampling, though it")
        print("still falls short of (does not overlap) the literature polarization-switching")
        print("window. This is a real, if modest, negative result for the bell hypothesis at this")
        print("specific location -- not an artifact of the fitting procedure.")
    else:
        print(f"S_peak std={S_peak_std:.3g} is intermediate: not tightly clustered, but not")
        print("scattered across the entire parameter range either.")
        print("\nVERDICT: Mixed/marginal robustness. Treat S_peak=1.4 as suggestive at best, not")
        print("confirmed -- resampling neither strongly supports nor strongly refutes it.")

    print(f"\nFor reference: bootstrap S_peak 90% interval = [{S_peak_p5:.4g}, {S_peak_p95:.4g}]")
    print(f"width 90% interval = [{width_p5:.4g}, {width_p95:.4g}]")
    print(f"kappa_max 90% interval = [{kappa_p5:.4g}, {kappa_p95:.4g}]")

    n_loo_failed = sum(1 for r in loo_results if r is None)
    n_loo_at_bound = sum(1 for r in loo_results if r is not None and any_at_bound(r))
    print(f"\nLeave-one-out: {N - n_loo_failed}/{N} converged, "
          f"{n_loo_at_bound}/{max(N - n_loo_failed, 1)} of those at a parameter bound -- "
          f"{'consistent with' if (n_loo_failed > 0 or n_loo_at_bound > 0) else 'not inconsistent with'} "
          f"the bootstrap picture above.")


if __name__ == "__main__":
    main()

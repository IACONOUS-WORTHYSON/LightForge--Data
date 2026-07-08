#!/usr/bin/env python3
"""
compare_light_field_models.py
==============================

Runs all three light-field functional-form calibrations against the SAME real
glycine NPLIN data (Javid et al. 2016, Cryst. Growth Des. 16, 4196-4202; see
GLYCINE_JAVID2016 in each model file for the full citation), overlays their
nucleation-rate ratio curves on one plot, and prints an honest plain-language
comparison of which (if any) predict a detectable effect in the real
experimental range S=1.4-1.7.

The three models compared (each a separate, independently rerunnable script --
see that script for its own full physics writeup and console diagnostics):

    1. nplin_cnt_model.py            -- ORIGINAL: light field subtracts from the
       nucleation barrier dG* (fixed coupling, independent of S).
    2. model_supersaturation_variant.py -- kappa_effective = kappa*(S-1): coupling
       hypothesized to grow with supersaturation (hypothesis, not physics).
    3. model_prefactor_variant.py    -- light field boosts the kinetic prefactor A
       instead of the barrier. IMPORTANT: this one is mathematically degenerate
       with A at the single fixed laser intensity in the real calibration data --
       its curve here is the upper limit the plausibility bound permits, NOT a
       fitted result (see that script's own output for the full explanation).

Run this after (or instead of) running the three scripts' own `calibrate-glycine`
mode individually -- each of those still works standalone and prints far more
detail about its own fit. This script's job is only the side-by-side comparison.

Usage:
    python compare_light_field_models.py
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def _load_module(filename, module_name):
    """Load a sibling script as an importable module without running its
    `if __name__ == "__main__"` block (that block is only entered when the file
    is executed directly, not when imported like this)."""
    spec = importlib.util.spec_from_file_location(module_name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    # Must register in sys.modules BEFORE exec_module: the @dataclass decorator
    # (used by MaterialParams etc.) looks up its own module via sys.modules to
    # resolve type hints, and fails with a confusing AttributeError if the
    # module isn't registered yet.
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _run_quietly(fn):
    """Run fn() with its stdout captured (each model's calibrate-glycine mode is
    very verbose -- that detail belongs to running that script directly, not to
    this comparison's console output) and return (result, captured_text)."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        result = fn()
    return result, buf.getvalue()


def main():
    print("Loading and running all three models' calibrate-glycine fits "
          "(each against the same real Javid et al. 2016 data)...")
    print("(Full diagnostic output from each model is suppressed here for readability --")
    print(" rerun any of the three scripts directly with 'calibrate-glycine' to see it.)\n")

    original = _load_module("nplin_cnt_model.py", "nplin_original")
    s_scaling = _load_module("model_supersaturation_variant.py", "nplin_sscaling")
    prefactor = _load_module("model_prefactor_variant.py", "nplin_prefactor")

    results = {}
    logs = {}
    for key, mod in [("original", original), ("s_scaling", s_scaling), ("prefactor", prefactor)]:
        print(f"  running {mod.__name__}...")
        result, log = _run_quietly(mod.run_calibrate_glycine)
        results[key] = result
        logs[key] = log

    print("\nAll three fits complete. Building overlay plot...\n")

    # ------------------------------------------------------------------
    # Overlay plot
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5.5))
    colors = {"original": "tab:purple", "s_scaling": "tab:blue", "prefactor": "tab:orange"}
    styles = {"original": "-", "s_scaling": "-", "prefactor": "--"}

    for key in ["original", "s_scaling", "prefactor"]:
        r = results[key]
        ax.semilogy(r["S_range"], r["ratio"], styles[key], color=colors[key], lw=2.2,
                    label=f"{r['label']} (kappa={r['kappa']:.3g}, R^2={r['r_squared']:.2f})")

    ax.axhline(1.0, color="k", ls=":", lw=1, label="no effect (ratio = 1)")
    ax.axvspan(1.4, 1.7, color="tab:green", alpha=0.12, label="target range S=1.4-1.7")
    ax.set_xlabel("Supersaturation, S")
    ax.set_ylabel("Nucleation rate ratio, J(with field) / J(no field)")
    ax.set_title("Three light-field functional forms, calibrated against the SAME real\n"
                 "glycine NPLIN data (Javid et al. 2016) -- NOT diamond/carbon/CO2\n"
                 "Prefactor curve = plausibility-bound UPPER LIMIT, not a fit (kappa\n"
                 "unidentifiable there -- see its own script's output)",
                 fontsize=9)
    ax.legend(fontsize=7.5, loc="upper right")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    out_path = HERE / "compare_light_field_models.png"
    fig.savefig(out_path, dpi=150)
    print(f"Saved comparison plot to {out_path.name}")

    # ------------------------------------------------------------------
    # Quantitative comparison within the real experimental range
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("QUANTITATIVE COMPARISON WITHIN S=1.4-1.7 (the real NPLIN experimental range)")
    print("=" * 78)
    print(f"{'Model':<45} {'ratio range in S=1.4-1.7':<28} {'max dev. from 1':>15}")
    summaries = {}
    for key in ["original", "s_scaling", "prefactor"]:
        r = results[key]
        mask = (r["S_range"] >= 1.4) & (r["S_range"] <= 1.7)
        r_lo, r_hi = r["ratio"][mask].min(), r["ratio"][mask].max()
        dev_pct = max(abs(r_lo - 1), abs(r_hi - 1)) * 100
        summaries[key] = (r_lo, r_hi, dev_pct)
        print(f"{r['label']:<45} {f'{r_lo:.3g} - {r_hi:.3g}':<28} {dev_pct:>14.3g}%")

    # ------------------------------------------------------------------
    # Plain-language summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("PLAIN-LANGUAGE SUMMARY")
    print("=" * 78)
    print("All three models were calibrated against the SAME 3 real glycine data points")
    print("(Javid et al. 2016), under the SAME literature-plausibility bound (rate ratio")
    print("<= 100x anywhere in S=1.05-2.0), fit with kappa properly search-bounded during")
    print("optimization (not clipped afterward -- see each script's own output).\n")

    orig_dev = summaries["original"][2]
    sscale_dev = summaries["s_scaling"][2]
    prefactor_dev = summaries["prefactor"][2]

    print(f"- ORIGINAL (fixed barrier subtraction): {orig_dev:.2g}% deviation in S=1.4-1.7.")
    if orig_dev < 1.0:
        print("  Negligible -- this mechanism's leverage is concentrated near the classical")
        print("  nucleation threshold (S~1.0-1.2) and has decayed away by S=1.4.")
    elif orig_dev < 5.0:
        print("  Small (1-5%) -- probably undetectable against realistic trial-to-trial noise.")
    else:
        print("  Meaningful (>=5%) within the real experimental range.")

    print(f"\n- SUPERSATURATION-SCALING (kappa*(S-1)): {sscale_dev:.2g}% deviation in S=1.4-1.7.")
    if sscale_dev < 1.0:
        print("  Also negligible -- even though this hypothesis makes the coupling grow with S,")
        print("  it is not enough to overcome how fast the underlying CNT barrier itself")
        print("  flattens out away from the threshold. The functional form matters less than")
        print("  the underlying (ln S)^-2 barrier shape, which both barrier-based models share.")
    elif sscale_dev < 5.0:
        print("  Small (1-5%) -- somewhat more persistent than the original fixed-coupling form,")
        print("  but still likely undetectable against realistic trial-to-trial noise.")
    else:
        print("  Meaningfully larger than the original fixed-coupling model in this range --")
        print("  the growing-with-S hypothesis does help the effect persist here.")

    print(f"\n- PRE-EXPONENTIAL-FACTOR (kappa boosts A): shows {prefactor_dev:.2g}% deviation,")
    print("  but THIS NUMBER IS NOT A FITTED RESULT. kappa is mathematically degenerate with A")
    print("  at the single fixed laser intensity this real dataset uses -- any kappa from 0 up")
    print("  to the plausibility bound fits equally well (confirmed in that script's own output:")
    print("  R^2 barely changes between kappa=0 and kappa=kappa_max). The curve shown is the")
    print("  most optimistic case the bound allows, not evidence the data prefer it over zero.")
    print("  This functional form cannot be tested with this dataset at all -- it would need")
    print("  nucleation data at multiple laser intensities, not just multiple S at one intensity.")

    print("\n--- BOTTOM LINE ---")
    testable_devs = {"original": orig_dev, "s_scaling": sscale_dev}
    best_key = max(testable_devs, key=testable_devs.get)
    best_dev = testable_devs[best_key]
    if best_dev < 1.0:
        print("Of the two functional forms this real dataset can actually test (original and")
        print("supersaturation-scaling), NEITHER predicts a detectable effect in the real")
        print("S=1.4-1.7 NPLIN experimental range once kappa is held to a literature-plausible")
        print("rate-ratio bound. Both concentrate their leverage very close to the classical")
        print("nucleation threshold (S~1.0-1.2), well below where real glycine NPLIN experiments")
        print("are run. The third form (prefactor) is simply untestable with this dataset.")
        print("Honest conclusion: this real glycine calibration gives NO support, under any of")
        print("the three functional forms tried, for a physically-plausible light-field effect")
        print("persisting at realistic NPLIN supersaturations. This does not disprove NPLIN --")
        print("it says a small, uniform free-energy-perturbation mechanism of the kind modeled")
        print("here is not what's needed to explain the size of the real effect in this dataset;")
        print("recall the source paper's own conclusion points at impurity/thermocavitation-")
        print("mediated nucleation instead, which none of these three variants model.")
    elif best_dev < 5.0:
        print(f"The best-performing testable model ({results[best_key]['label']}) gives only a")
        print(f"small ({best_dev:.2g}%) effect in S=1.4-1.7 -- likely undetectable against realistic")
        print("trial-to-trial nucleation-probability noise without a very large repeat-trial count.")
        print("This is weak, inconclusive support at best, not a confirmed effect.")
    else:
        print(f"The best-performing testable model ({results[best_key]['label']}) gives a")
        print(f"meaningful ({best_dev:.2g}%) effect within S=1.4-1.7 under the plausibility bound.")
        print("This is the most encouraging result among the three, though still only as strong")
        print("as a 3-point calibration with unmeasured placeholder gamma/A can support.")


if __name__ == "__main__":
    main()

"""
Diamond CVD Integrated Fabrication-Time Model
==============================================

Estimates total fabrication time (nucleation/incubation + steady-state
growth) for a target diamond part grown by microwave-plasma chemical
vapor deposition (MPCVD), optionally combined with two INDEPENDENTLY
PUBLISHED enhancement techniques:

    1. UV-laser-assisted growth (Lu et al. 2017)          -> growth-rate boost
    2. Femtosecond-laser nucleation seeding (review lit.)  -> shorter incubation
    3. Femtosecond-laser defect correction (2020 paper)    -> quality metric only,
       tracked separately, never folded into the time numbers.
    4. Masked selective-area growth (standard fab technique) -> defines WHICH
       volume/thickness must be grown, not the rate itself.

This script models a PROPOSED NOVEL COMBINATION of these techniques.
No paper has built or tested this combination as an integrated system.
Every numeric assumption below is individually commented as either a
direct literature citation or an explicitly-labeled estimate.

IMPORTANT SCOPE NOTE:
This model does NOT include, assume, or reference any "nucleation-bias"
/ kappa-style light-field effect. That was a separate hypothesis, already
tested independently, and returned a null (negligible) result. It is
intentionally excluded here and is unrelated to the techniques modeled
below (all of which are confirmed, published effects).

Literature referenced:
  [A] General MPCVD growth-rate baseline (1-10 um/hr for high-quality
      single-crystal diamond): widely-reported consensus range in CVD
      diamond literature, e.g. May, P.W., "Diamond thin films: a
      21st-century material," Phil. Trans. R. Soc. Lond. A 358 (2000);
      Butler, J.E. & Sumant, A.V., "The CVD of Nanodiamond Materials,"
      Chem. Vap. Deposition 14 (2008). Not a single-paper number.
  [B] Lu et al., "Ultraviolet laser photolysis of hydrocarbons for
      nondiamond carbon suppression in chemical vapor deposition of
      diamond films," Light: Science & Applications (2017), PMC6060054.
  [C] "Improvement of crystallization in CVD diamond coating induced by
      femtosecond laser irradiation," Diamond and Related Materials (2020).
  [D] Femtosecond-laser-induced nucleation seeding review/study on
      heteroepitaxial diamond growth (higher, more uniform nucleation
      density vs. mechanical-grinding seeding; exact time savings not
      quantified in source -> modeled as a conservative estimate).
  [E] Masked selective-area deposition: standard materials-fabrication
      technique, not diamond-specific literature.

Rerun this script directly: `python diamond_cvd_fabrication_model.py`
"""

import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# LITERATURE-SOURCED CONSTANTS AND EXPLICITLY-LABELED ESTIMATES
# ---------------------------------------------------------------------------

# [A] Baseline steady-state growth-rate range for high-quality single-crystal
# MPCVD diamond. Widely-reported consensus baseline, general CVD diamond
# literature (see module docstring, ref [A]). CONFIRMED RANGE.
BASELINE_RATE_MIN_UM_HR = 1.0
BASELINE_RATE_MAX_UM_HR = 10.0

# [B] UV-assist growth-rate multiplier. Lu et al. (2017), PMC6060054,
# report UV laser irradiation of the CVD flame DOUBLED diamond growth
# rate relative to standard CVD. CONFIRMED, DIRECTLY CITED NUMBER.
UV_ASSIST_GROWTH_MULTIPLIER = 2.0

# [B] UV-assist film-quality improvement, informational only -- NOT used
# in the time calculation. Lu et al. (2017) report ~4.2% quality
# improvement under UV assist. CONFIRMED, DIRECTLY CITED NUMBER.
UV_ASSIST_QUALITY_IMPROVEMENT_PCT = 4.2

# [B] Lu et al. (2017) also report UV assist "shortened nucleation time"
# qualitatively, but give no percentage/quantity. To avoid injecting an
# unsourced number into the model, this effect is NOT quantified here --
# it is reported as a qualitative note only in the summary text below.
UV_NUCLEATION_EFFECT_QUANTIFIED = False  # intentionally not modeled numerically

# [C] Femtosecond-laser defect-correction crystallinity improvement.
# "Improvement of crystallization in CVD diamond coating induced by
# femtosecond laser irradiation," Diamond and Related Materials (2020):
# ~20% crystallinity improvement at low fluence via photo-induced
# amorphous-to-diamond phase conversion. CONFIRMED, DIRECTLY CITED NUMBER.
# Modeled strictly as a QUALITY metric, reported separately -- it is never
# folded into the growth-rate or fabrication-time numbers, per the source
# describing it as a concurrent/post-growth correction step, not a rate
# multiplier.
FS_DEFECT_CORRECTION_CRYSTALLINITY_IMPROVEMENT_PCT = 20.0

# ESTIMATED, NOT DIRECTLY SOURCED. Representative order-of-magnitude
# placeholder for the nucleation/incubation delay before steady-state
# growth begins in standard MPCVD with conventional (e.g. mechanical
# grinding/scratch) seeding. Commonly described in the literature as
# ranging from tens of minutes to a few hours depending on substrate
# preparation and seeding density, but no single hard number is asserted
# here. User-adjustable; treat as a stand-in default, not a citation.
DEFAULT_NUCLEATION_HOURS = 1.0

# [D] ESTIMATED, NOT DIRECTLY SOURCED -- CONSERVATIVE PLACEHOLDER.
# Femtosecond-laser-induced nucleation seeding is reported to produce
# higher nucleation density and more uniform distribution than mechanical-
# grinding seeding for heteroepitaxial diamond growth (ref [D]), which
# qualitatively implies a shorter incubation period before steady-state
# growth. The source material does NOT quantify the time savings, so this
# is a conservative 50% reduction assumption used purely for modeling
# purposes -- flag clearly as an estimate whenever quoting this figure.
FS_SEEDING_NUCLEATION_REDUCTION_FRACTION = 0.5

# [E] Masked selective-area growth is a standard fabrication technique:
# growth only occurs where the mask exposes the substrate, and all
# exposed regions grow SIMULTANEOUSLY under the same plasma exposure.
# ENGINEERING ASSUMPTION (not a literature citation, just deposition
# physics): total fabrication time is governed by the TALLEST required
# thickness in the shape's thickness map, not the sum of all volumes --
# thin and thick regions grow in parallel, so extra volume elsewhere
# doesn't add extra time, it only adds extra deposited material.
MASK_TIME_GOVERNED_BY = "max_thickness"  # documents the assumption above


# ---------------------------------------------------------------------------
# CORE MODEL
# ---------------------------------------------------------------------------

def steady_state_growth_hours(thickness_um, effective_rate_um_hr):
    """Hours of steady-state growth to reach a given thickness at a given rate."""
    return thickness_um / effective_rate_um_hr


def nucleation_phase_hours(base_nucleation_hours=DEFAULT_NUCLEATION_HOURS,
                            fs_seeding=False):
    """
    Nucleation/incubation time before steady-state growth begins.

    fs_seeding=True applies the conservative, explicitly-estimated
    FS_SEEDING_NUCLEATION_REDUCTION_FRACTION (see comment above, source [D]).
    """
    if fs_seeding:
        return base_nucleation_hours * (1.0 - FS_SEEDING_NUCLEATION_REDUCTION_FRACTION)
    return base_nucleation_hours


def effective_growth_rate(baseline_rate_um_hr, uv_assist=False):
    """
    Applies the CONFIRMED UV_ASSIST_GROWTH_MULTIPLIER (source [B]) to the
    chosen baseline rate when uv_assist is enabled.
    """
    if not (BASELINE_RATE_MIN_UM_HR <= baseline_rate_um_hr <= BASELINE_RATE_MAX_UM_HR):
        raise ValueError(
            f"baseline_rate_um_hr must be within the literature range "
            f"[{BASELINE_RATE_MIN_UM_HR}, {BASELINE_RATE_MAX_UM_HR}] um/hr (source [A])."
        )
    multiplier = UV_ASSIST_GROWTH_MULTIPLIER if uv_assist else 1.0
    return baseline_rate_um_hr * multiplier


def required_deposition_volume_mm3(thickness_map):
    """
    Informational only: total material volume that must be deposited
    through the mask openings, given a thickness map
    [{"thickness_um": ..., "area_mm2": ...}, ...] (source [E]).
    Does NOT affect fabrication time in this model (see MASK_TIME_GOVERNED_BY).
    """
    total_mm3 = 0.0
    for region in thickness_map:
        thickness_mm = region["thickness_um"] * 1e-3
        total_mm3 += thickness_mm * region["area_mm2"]
    return total_mm3


def max_required_thickness_um(thickness_map):
    """The tallest required feature governs total growth time (source [E] + parallel-growth assumption)."""
    return max(region["thickness_um"] for region in thickness_map)


def estimate_fabrication_time(thickness_map,
                               baseline_rate_um_hr,
                               uv_assist=False,
                               fs_seeding=False,
                               base_nucleation_hours=DEFAULT_NUCLEATION_HOURS):
    """
    Full fabrication-time estimate for a masked selective-area target shape.

    thickness_map: list of {"thickness_um": float, "area_mm2": float}
                    (for a first-version single-thickness part, pass a
                    single-entry list.)
    Returns a dict with nucleation, growth, and total hours, plus
    informational effective rate and deposition volume.
    """
    governing_thickness_um = max_required_thickness_um(thickness_map)
    rate = effective_growth_rate(baseline_rate_um_hr, uv_assist=uv_assist)

    nuc_h = nucleation_phase_hours(base_nucleation_hours, fs_seeding=fs_seeding)
    growth_h = steady_state_growth_hours(governing_thickness_um, rate)
    total_h = nuc_h + growth_h

    return {
        "governing_thickness_um": governing_thickness_um,
        "effective_rate_um_hr": rate,
        "nucleation_hours": nuc_h,
        "growth_hours": growth_h,
        "total_hours": total_h,
        "deposition_volume_mm3": required_deposition_volume_mm3(thickness_map),
    }


# ---------------------------------------------------------------------------
# DEMO / CONFIG -- edit these to explore other scenarios
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    BASELINE_RATE_UM_HR = 5.0          # pick any value within [1, 10] um/hr, source [A]
    EXAMPLE_AREA_MM2 = 100.0           # informational, affects volume not time
    EXAMPLE_THICKNESSES_UM = [50, 200, 500]

    scenarios = {
        "Standard CVD only":        dict(uv_assist=False, fs_seeding=False),
        "UV-assisted":              dict(uv_assist=True,  fs_seeding=False),
        "UV-assisted + FS seeding": dict(uv_assist=True,  fs_seeding=True),
    }

    print("=" * 78)
    print("DIAMOND CVD INTEGRATED FABRICATION-TIME MODEL")
    print("=" * 78)
    print(f"Baseline growth rate used: {BASELINE_RATE_UM_HR} um/hr "
          f"(within literature range {BASELINE_RATE_MIN_UM_HR}-{BASELINE_RATE_MAX_UM_HR} um/hr, source [A])")
    print(f"Default nucleation/incubation baseline: {DEFAULT_NUCLEATION_HOURS} hr "
          f"(ESTIMATED, not directly sourced)")
    print()

    results = {name: [] for name in scenarios}

    for thickness_um in EXAMPLE_THICKNESSES_UM:
        thickness_map = [{"thickness_um": thickness_um, "area_mm2": EXAMPLE_AREA_MM2}]
        print(f"--- Target thickness: {thickness_um} um  (area {EXAMPLE_AREA_MM2} mm^2) ---")
        for name, opts in scenarios.items():
            r = estimate_fabrication_time(
                thickness_map,
                baseline_rate_um_hr=BASELINE_RATE_UM_HR,
                **opts,
            )
            results[name].append(r["total_hours"])
            print(f"  {name:28s}: nucleation {r['nucleation_hours']:5.2f} hr "
                  f"+ growth {r['growth_hours']:6.2f} hr "
                  f"= total {r['total_hours']:6.2f} hr "
                  f"(effective rate {r['effective_rate_um_hr']:.1f} um/hr)")
        print()

    # -----------------------------------------------------------------------
    # CHART: standard vs UV-assisted (vs UV+FS-seeding) fabrication time
    # -----------------------------------------------------------------------
    x = range(len(EXAMPLE_THICKNESSES_UM))
    bar_width = 0.25

    fig, ax = plt.subplots(figsize=(8, 5))
    for i, (name, totals) in enumerate(results.items()):
        offset = (i - 1) * bar_width
        positions = [xi + offset for xi in x]
        ax.bar(positions, totals, width=bar_width, label=name)

    ax.set_xticks(list(x))
    ax.set_xticklabels([f"{t} um" for t in EXAMPLE_THICKNESSES_UM])
    ax.set_xlabel("Target thickness")
    ax.set_ylabel("Estimated fabrication time (hours)")
    ax.set_title("Diamond CVD Fabrication Time: Standard vs. UV-Assisted\n"
                  f"(baseline rate {BASELINE_RATE_UM_HR} um/hr, area {EXAMPLE_AREA_MM2} mm^2)")
    ax.legend()
    fig.tight_layout()

    output_path = "diamond_cvd_fabrication_time_comparison.png"
    fig.savefig(output_path, dpi=150)
    print(f"Chart saved to: {output_path}")

    # -----------------------------------------------------------------------
    # PLAIN-LANGUAGE SUMMARY
    # -----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print(
        "This model projects diamond CVD fabrication time using CONFIRMED,\n"
        "PUBLISHED growth-rate and process improvements only:\n"
        "  - Baseline MPCVD growth rate: 1-10 um/hr (general literature range).\n"
        "  - UV-assisted growth: 2x growth-rate multiplier, directly cited from\n"
        "    Lu et al. (2017), Light: Science & Applications, PMC6060054.\n"
        "  - Femtosecond-laser nucleation seeding: incubation-time reduction is\n"
        "    an EXPLICITLY LABELED, CONSERVATIVE ESTIMATE (50%), since the\n"
        "    source material reports higher/more-uniform nucleation density but\n"
        "    does not quantify the resulting time savings.\n"
        "  - Femtosecond-laser defect correction (~20% crystallinity improvement,\n"
        "    Diamond and Related Materials, 2020) is tracked as a QUALITY metric\n"
        "    only and is deliberately NOT folded into any growth-rate or\n"
        "    fabrication-time number.\n"
        "  - Masked selective-area growth defines WHICH volume must be grown; it\n"
        "    does not change the rate itself.\n"
        "\n"
        "This model does NOT include or assume any nucleation-bias / kappa-style\n"
        "light-field effect. That mechanism was separately tested and found\n"
        "negligible at realistic conditions, and is intentionally excluded here.\n"
        "\n"
        "This model represents a PROPOSED NOVEL COMBINATION of individually-\n"
        "published techniques. No paper or prototype has built or tested this\n"
        "combination as an integrated system -- these are projected estimates,\n"
        "not measured results."
    )

#!/usr/bin/env python3
"""
nplin_cnt_model.py
===================

Classical Nucleation Theory (CNT) simulator with an optional tunable
"light-field" term motivated by the non-photochemical laser-induced
nucleation (NPLIN) literature, a Monte Carlo stochastic layer, and a
calibration mode for fitting against real nucleation-probability-vs-fluence
data.

THIS IS A PLAUSIBILITY / PARAMETER-EXPLORATION TOOL, NOT A VALIDATED
PREDICTIVE MODEL. NPLIN's physical mechanism is still actively debated in
the literature (see "Physics background" below) -- use this script to sanity
check orders of magnitude and to plan experiments, then calibrate against
your own or published data before trusting any quantitative prediction.

-----------------------------------------------------------------------------
PHYSICS BACKGROUND AND REFERENCES
-----------------------------------------------------------------------------

Classical nucleation theory (the part with no light field):

    Kashchiev, D. "Nucleation: Basic Theory with Applications."
    Butterworth-Heinemann, 2000.  (source of the standard CNT formulas used
    below: barrier height, critical radius, Arrhenius-form rate law.)

    Mullin, J. W. "Crystallization," 4th ed. Butterworth-Heinemann, 2001.
    (general source for typical orders of magnitude of the CNT kinetic
    prefactor A; this prefactor is notoriously uncertain -- prefer fitting
    it to your own induction-time or calibration data over trusting any
    textbook default.)

NPLIN phenomenology and proposed mechanisms (the light-field term below is
a simplified, tunable implementation of ONE of at least three mechanisms
that have been proposed in this literature -- see caveats further down):

    Garetz, B. A.; Aber, J. E.; Goddard, N. L.; Young, R. G.; Myerson, A. S.
    "Nonphotochemical, Polarization-Dependent, Laser-Induced Nucleation in
    Supersaturated Aqueous Urea Solutions." Phys. Rev. Lett. 1996, 77, 3475.
    (original NPLIN report; proposed the optical-Kerr / anisotropic-molecule
    alignment mechanism.)

    Garetz, B. A.; Matic, J.; Myerson, A. S. "Polarization Switching of
    Crystal Structure in the Nonphotochemical Light-Induced Nucleation of
    Supersaturated Aqueous Glycine Solutions." Phys. Rev. Lett. 2002, 89,
    175501. (glycine polymorph switching with polarization -- relevant if
    you are using glycine as your model compound.)

    Knott, B. C.; Doherty, M. F.; Peters, B. "A simulation test of the
    optical Kerr mechanism for laser-induced nucleation." J. Chem. Phys.
    2011, 134, 154501. (Monte Carlo lattice-gas test showing the Kerr
    /orientational-alignment mechanism needs field strengths far above
    experimental values to explain observed nucleation enhancement --
    motivated the search for alternative mechanisms.)

    Alexander, A. J.; Camp, P. J. "Non-photochemical laser-induced
    nucleation." J. Chem. Phys. 2019, 150, 040901 (review/perspective).
    Proposes and reviews the "Dielectric Polarization" (DP) model: even
    ISOTROPIC clusters (e.g. cubic KCl, which has no anisotropic molecules
    to align) can be stabilized by an optical field if the cluster's
    (electronic, i.e. optical-frequency) dielectric constant differs from
    that of the surrounding solvent -- a dielectric sphere in an external
    field lowers its energy, which lowers the free energy of pre-critical
    clusters and can shrink the critical nucleus size. THIS is the
    mechanism the light-field term in this script implements, because it
    is the one explicitly described as an "isotropic electronic
    polarization" effect in the request this script was built for.

    Garetz, B. A.; Hartman, R. L. "25+ Years of Research on Nonphotochemical
    Laser-Induced Nucleation (NPLIN)." Cryst. Growth Des. 2025, 25, 2756.
    (recent perspective summarizing all three proposed mechanism families:
    optical Kerr / orientational, dielectric polarization / isotropic
    electronic, and impurity-mediated thermocavitation.)

    Electrostriction / dielectric-sphere-in-a-field energy expression:
    Landau, L. D.; Lifshitz, E. M. "Electrodynamics of Continuous Media,"
    2nd ed., Pergamon, 1984 (standard electromagnetism result for the energy
    of a dielectric sphere embedded in a different dielectric medium under
    an applied field; used here in its low local-field / Clausius-Mossotti
    form).

CAVEAT YOU SHOULD NOT SKIP: a growing body of more recent NPLIN work (e.g.
the cavitation/thermocavitation studies discussed in Garetz & Hartman 2025,
and dedicated impurity-nanoparticle studies such as the CrystEngComm 2024
sodium-acetate cavitation paper, DOI 10.1039/D4CE00487F) argues that in many
published NPLIN experiments the dominant effect is NOT direct field-cluster
coupling at all, but laser heating of solid impurity nanoparticles causing
thermocavitation (a transient vapor bubble) that nucleates the crystal by a
completely different, impurity-mediated pathway. That mechanism is NOT
modeled here (it depends on impurity loading, filtration, particle
absorption cross-sections -- not on the clean-solution electronic-structure
picture this script implements). If your diamond-fabrication concept
depends specifically on DIRECT optical-field/cluster coupling (as opposed to
"laser pulses happen to nucleate stuff, by whatever mechanism"), treat any
encouraging fold-change this script predicts as an upper-bound / best-case
plausibility check for that specific hypothesis, not as evidence that the
coupling pathway is the one actually responsible for any nucleation you see
experimentally. Filtering your solution and re-running is the standard way
experimentalists in this field distinguish the two.

-----------------------------------------------------------------------------
CNT MODEL USED
-----------------------------------------------------------------------------

For a spherical nucleus of radius r growing from solution at supersaturation
S = c / c_eq:

    dG(r) = 4*pi*r^2*gamma - (4/3)*pi*r^3*(kT/v_m)*ln(S)

    r*    = 2*gamma*v_m / (kT*ln S)                      (critical radius)
    dG*   = 16*pi*gamma^3*v_m^2 / (3*(kT)^2*(ln S)^2)     (barrier height)
    J     = A * exp(-dG*/kT)                              (nucleation rate,
                                                            events per m^3
                                                            per s)

gamma = solid-liquid interfacial free energy [J/m^2]
v_m   = volume per molecule in the solid [m^3]
A     = kinetic prefactor [m^-3 s^-1] (attachment rate x monomer density x
        Zeldovich factor, lumped; treat as an empirical/fit parameter)

-----------------------------------------------------------------------------
LIGHT-FIELD TERM (DP / isotropic electronic polarization mechanism)
-----------------------------------------------------------------------------

A dielectric sphere of (optical-frequency) permittivity eps_c, volume V,
embedded in a medium of permittivity eps_s, in a field of peak amplitude E0,
has field-induced energy (Clausius-Mossotti / low-contrast approximation):

    u_field = (1/2) * eps0 * K(eps_c, eps_s) * E0^2      [J/m^3, i.e. energy
                                                           per unit cluster
                                                           volume]
    K(eps_c, eps_s) = 3*eps_s*(eps_c - eps_s) / (eps_c + 2*eps_s)

The FREE-ENERGY-OF-FORMATION contribution of this interaction (the quantity
that adds to the bulk driving-force term below) is the negative of this
interaction energy density, i.e. it is POSITIVE and increases the effective
driving force when eps_c > eps_s -- e.g. when the pre-nucleating cluster is
optically denser than the surrounding solution, which is the qualitative
picture in Alexander & Camp's DP model. (Equivalently: the field lowers the
cluster's own energy, which is a lowering of the free-energy cost of forming
it -- a bigger driving force, not a smaller one.) We treat eps as the
optical-frequency (electronic) permittivity,
eps ~= n^2, since only the fast electronic polarizability follows an
optical-frequency E-field -- this is what makes the mechanism apply even to
isotropic/non-birefringent crystals, unlike the orientational optical-Kerr
picture.

We fold ALL of the remaining physical uncertainty (how much of this nominal
per-cluster energy actually acts on a real pre-critical cluster; local-field
enhancement; deviation from a perfect sphere; whether subcritical clusters
even have bulk-crystal-like permittivity, etc.) into one dimensionless,
user-tunable "coupling efficiency" kappa. kappa = 0 recovers ordinary CNT.
kappa = 1 means "take the nominal dielectric-sphere energy at face value."
kappa is exactly the parameter you should scan to test different hypotheses
about how strongly light couples to nucleation in your system, and it is
also the parameter the calibrate mode fits against real data.

The field energy adds directly to the bulk (volume) free-energy term,
because both scale with cluster volume:

    ln(S_eff) = ln(S) + kappa * u_field * v_m / kT

and dG*, r*, J are then evaluated at S_eff instead of S. u_field > 0 (for
eps_c > eps_s) increases ln(S_eff), which shrinks r* and dG* -- i.e. exactly
the "reduction in free energy of pre-nucleating clusters ... reducing
critical nucleus size" effect described in the NPLIN literature above.

-----------------------------------------------------------------------------
STOCHASTIC / MONTE CARLO LAYER
-----------------------------------------------------------------------------

CNT gives a RATE, not a certainty. For a solution volume V observed for a
time t, the probability of at least one nucleation event (Poisson process
with rate J*V) is:

    P = 1 - exp(-J*V*t)

For n_pulses repeated, independent pulses each of duration tau:

    P_total = 1 - (1 - P_single_pulse)^n_pulses

Real repeated trials of nominally identical experiments are then modeled as
independent Bernoulli(P) draws -- this is the same logic used to report
"nucleation probability" in real NPLIN papers (fraction of vials/drops that
nucleated out of N attempts), and is why NPLIN data is probabilistic across
repeated trials rather than deterministic.

-----------------------------------------------------------------------------
USAGE
-----------------------------------------------------------------------------

Requires: numpy, scipy, matplotlib (no other dependencies).

    python nplin_cnt_model.py demo
        Fully self-contained illustrative walkthrough using CLEARLY LABELED
        placeholder parameters and SYNTHETIC (not real) calibration data.
        Run this first to see the whole pipeline work end to end.

    python nplin_cnt_model.py explore
        Interactively prompts you for your system's real parameters, then
        plots nucleation rate vs. supersaturation with and without the
        light-field term, and prints the predicted fold-change and
        critical-radius reduction.

    python nplin_cnt_model.py calibrate --data yourdata.csv
        Loads real (fluence, nucleation probability[, n_trials]) data you
        have digitized from a published paper or measured yourself, and
        fits the coupling efficiency kappa (and optionally the CNT
        prefactor A) to match it, at a supersaturation you specify.
        CSV format, header required:
            fluence_J_cm2,probability[,n_trials]

    python nplin_cnt_model.py montecarlo
        Runs repeated stochastic trials (binomial sampling) across a range
        of supersaturations at your specified laser settings, and plots the
        simulated observed fraction (with Wilson confidence intervals)
        against the analytic probability curve, with and without the light
        field.

    python nplin_cnt_model.py litcheck
        Sanity-checks the Monte Carlo/probability layer against REAL published
        repeat-trial nucleation counts (KCl, KBr -- see citations in the code).
        These are not fluence sweeps and cannot calibrate kappa; this only
        validates the stochastic layer's statistics against real data.

Any of the three parameter-driven modes accept --config a_file.json to load
previously-entered parameters (skips the interactive prompts) and
--save-config a_file.json to write out whatever you just entered, so you can
rerun with tweaked inputs by hand-editing the JSON instead of re-typing
everything. This is meant to be reused, not run once.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
from scipy import constants as sc
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# =============================================================================
# Physical constants (SI)
# =============================================================================
K_B = sc.k            # Boltzmann constant, J/K
EPS0 = sc.epsilon_0    # vacuum permittivity, F/m
C_LIGHT = sc.c         # speed of light in vacuum, m/s

# Unit-conversion helpers (script works internally in SI; prompts accept the
# units people actually report things in)
MM_TO_M = 1e-3
NM_TO_M = 1e-9
NS_TO_S = 1e-9
UL_TO_M3 = 1e-9
MJ_TO_J = 1e-3
JCM2_TO_JM2 = 1e4


# =============================================================================
# Parameter containers
# =============================================================================
@dataclass
class MaterialParams:
    name: str
    gamma: float              # solid-liquid interfacial free energy, J/m^2
    v_m: float                # molecular volume in the solid, m^3/molecule
    T: float                  # temperature, K
    A_prefactor: float        # CNT kinetic prefactor, 1/(m^3 s)
    n_solvent: float          # refractive index of the mother liquor at the laser wavelength
    n_solid: float            # refractive index of the crystal/cluster at the laser wavelength


@dataclass
class LightFieldParams:
    wavelength_nm: float
    pulse_energy_J: float
    spot_diameter_m: float
    pulse_duration_s: float
    n_pulses: int
    coupling_efficiency: float   # kappa, dimensionless, tunable
    pulse_shape_factor: float = 1.0  # peak/average intensity correction; 1.0 = flat-top approx


@dataclass
class SampleGeometry:
    volume_m3: float          # illuminated / observed solution volume
    exposure_time_s: float    # holding/observation time per trial (post-pulse), used by nucleation_probability


# =============================================================================
# Core CNT model (no light field)
# =============================================================================
def delta_g_star(S, T, gamma, v_m):
    """Classical nucleation barrier height, J. See module docstring for the formula."""
    lnS = np.log(S)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 16 * np.pi * gamma**3 * v_m**2 / (3 * (K_B * T) ** 2 * lnS**2)


def critical_radius(S, T, gamma, v_m):
    """Critical nucleus radius, m."""
    lnS = np.log(S)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 2 * gamma * v_m / (K_B * T * lnS)


def nucleation_rate(S, T, gamma, v_m, A):
    """CNT nucleation rate, events m^-3 s^-1, with no light field."""
    dG = delta_g_star(S, T, gamma, v_m)
    return A * np.exp(-dG / (K_B * T))


# =============================================================================
# Light-field (DP / isotropic electronic polarization) term
# =============================================================================
def fluence_from_pulse(pulse_energy_J, spot_diameter_m):
    """Laser fluence (energy density), J/m^2, for a flat-top circular spot."""
    area = np.pi * (spot_diameter_m / 2) ** 2
    return pulse_energy_J / area


def peak_intensity_from_fluence(H_J_m2, pulse_duration_s, pulse_shape_factor=1.0):
    """Peak intensity, W/m^2, from fluence and pulse duration (flat-top approx by default)."""
    return pulse_shape_factor * H_J_m2 / pulse_duration_s


def dielectric_sphere_factor(eps_cluster, eps_solvent):
    """Clausius-Mossotti factor K for a dielectric sphere in a different dielectric medium.
    K > 0 (stabilizing) when eps_cluster > eps_solvent."""
    return 3 * eps_solvent * (eps_cluster - eps_solvent) / (eps_cluster + 2 * eps_solvent)


def field_energy_density(I_peak, n_solvent, n_solid, kappa):
    """Field-induced contribution to the bulk driving-force term, J/m^3 (this is
    -1x the dielectric-sphere interaction energy density; see module docstring
    for the distinction). Positive -- i.e. promotes nucleation, shrinks the
    barrier -- when n_solid > n_solvent. `kappa` is the tunable coupling-
    efficiency hypothesis parameter (see module docstring)."""
    eps_s = n_solvent**2
    eps_c = n_solid**2
    K = dielectric_sphere_factor(eps_c, eps_s)
    E0_sq = 2 * I_peak / (C_LIGHT * EPS0 * n_solvent)
    return 0.5 * EPS0 * K * E0_sq * kappa


def effective_ln_S(S, u_field, v_m, T):
    """Supersaturation term as modified by the light field's contribution to the
    bulk (volume) free-energy term."""
    return np.log(S) + u_field * v_m / (K_B * T)


def delta_g_star_field(S, T, gamma, v_m, u_field):
    kT = K_B * T
    lnS_eff = effective_ln_S(S, u_field, v_m, T)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 16 * np.pi * gamma**3 * v_m**2 / (3 * kT**2 * lnS_eff**2)


def critical_radius_field(S, T, gamma, v_m, u_field):
    lnS_eff = effective_ln_S(S, u_field, v_m, T)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 2 * gamma * v_m / (K_B * T * lnS_eff)


def nucleation_rate_field(S, T, gamma, v_m, A, u_field):
    """CNT nucleation rate including the light-field term."""
    dG = delta_g_star_field(S, T, gamma, v_m, u_field)
    return A * np.exp(-dG / (K_B * T))


def u_field_from_params(material: MaterialParams, light: LightFieldParams) -> float:
    """Convenience: compute u_field (J/m^3) directly from a pulse-energy/spot/duration
    description of the laser, as you'd specify from your experimental setup."""
    H = fluence_from_pulse(light.pulse_energy_J, light.spot_diameter_m)
    I_peak = peak_intensity_from_fluence(H, light.pulse_duration_s, light.pulse_shape_factor)
    return field_energy_density(I_peak, material.n_solvent, material.n_solid, light.coupling_efficiency)


# =============================================================================
# Stochastic / Monte Carlo layer
# =============================================================================
def nucleation_probability(J, volume_m3, time_s):
    """Probability of >=1 nucleation event (Poisson process), for a single exposure window."""
    return 1.0 - np.exp(-J * volume_m3 * time_s)


def probability_multi_pulse(J_per_pulse, volume_m3, pulse_duration_s, n_pulses):
    """Probability of >=1 nucleation event across n_pulses independent pulses,
    each contributing a nucleation rate J_per_pulse active for pulse_duration_s."""
    p_single = nucleation_probability(J_per_pulse, volume_m3, pulse_duration_s)
    return 1.0 - (1.0 - p_single) ** n_pulses


def monte_carlo_trials(P, n_trials, rng=None):
    """Simulate n_trials independent Bernoulli(P) nucleation attempts. Returns the
    array of 0/1 outcomes -- this is the "did it nucleate this time" stochastic layer."""
    rng = rng if rng is not None else np.random.default_rng()
    return rng.binomial(1, np.clip(P, 0.0, 1.0), size=n_trials)


def wilson_ci(k, n, z=1.96):
    """Wilson score confidence interval for a binomial proportion k/n (no extra deps)."""
    if n == 0:
        return (np.nan, np.nan)
    phat = k / n
    denom = 1 + z**2 / n
    center = (phat + z**2 / (2 * n)) / denom
    half = z * np.sqrt(phat * (1 - phat) / n + z**2 / (4 * n**2)) / denom
    return (center - half, center + half)


# =============================================================================
# Real literature trial-count data (repeat-trial validation of the stochastic layer)
# =============================================================================
# NOTE ON PROVENANCE: these numbers were extracted from the papers via an automated
# fetch-and-summarize tool, not read directly by a human off the original PDF/HTML.
# They are cited so you can independently verify them before relying on them for
# anything beyond a rough sanity check of the Monte Carlo / probability layer below.
# They are NOT a fluence sweep, so they cannot be used to fit the light-field
# coupling kappa -- see fit_calibration()/calibrate mode for that, which needs a
# real fluence-dependent probability dataset you supply.
#
# KCl repeatability data -- same nominal laser intensity and exposed volume, repeated
# in separate sessions (source paper reports this explicitly as a repeatability check):
#   "Effect of Laser-Exposed Volume and Irradiation Position on Nonphotochemical
#   Laser-Induced Nucleation of Potassium Chloride Solutions" (PMC10626568, open
#   access), fixed laser intensity 10-23 MW/cm^2.
KCL_REPEATABILITY_PMC10626568 = {
    "citation": ("Effect of Laser-Exposed Volume and Irradiation Position on "
                 "Nonphotochemical Laser-Induced Nucleation of Potassium Chloride "
                 "Solutions, PMC10626568 (open access), fixed ~10-23 MW/cm^2 intensity"),
    "note": "Three repeat sessions at the same nominal condition -- suitable for testing whether the observed spread in nucleation probability is consistent with ordinary binomial sampling noise around one shared probability.",
    "trials": [("session 1", 39, 95), ("session 2", 56, 99), ("session 3", 48, 97)],
}

# KBr data -- DIFFERENT solution-additive conditions (not repeats of one condition), at
# fixed 1.0 mJ/pulse, 60.7 MW/cm^2 peak intensity:
#   Li, S.; Xie, X.; Liu, Y. "Effect of acidic polymers on the morphology of
#   non-photochemical laser-induced nucleation of potassium bromide."
#   Sci. Rep. 2024, 14 (PMC10997761, open access).
KBR_ADDITIVES_PMC10997761 = {
    "citation": ("Li, Xie & Liu, \"Effect of acidic polymers on the morphology of "
                 "non-photochemical laser-induced nucleation of potassium bromide,\" "
                 "Sci. Rep. 2024, PMC10997761 (open access), fixed 1.0 mJ/pulse, "
                 "60.7 MW/cm^2 intensity"),
    "note": "Three DIFFERENT additive conditions, not repeats -- the paper's own point is that these probabilities differ. Do not pool these into a single-p consistency test.",
    "trials": [("no additive", 34, 35), ("19 wt% PAA", 20, 35), ("19 wt% PMA", 27, 35)],
}


def bootstrap_spread_pvalue(trials, n_boot=20000, rng=None):
    """Parametric bootstrap test: are these (k, n) repeat trials consistent with all being
    binomial draws around one shared pooled probability, or is the observed spread bigger
    than ordinary binomial sampling noise would produce? Reuses monte_carlo_trials(), the
    same stochastic layer used everywhere else in this script -- this is the same logic,
    just pointed at real repeat-trial counts instead of a synthetic example.

    Returns (pooled_p, observed_spread, p_value). A small p_value means the real spread is
    larger than binomial noise alone would typically produce (session-to-session physical
    variability, e.g. impurity content); a large p_value means pure binomial noise already
    explains it, with no need to invoke anything else.
    """
    rng = rng if rng is not None else np.random.default_rng()
    ks = np.array([k for _, k, n in trials], dtype=float)
    ns = np.array([n for _, k, n in trials], dtype=float)
    pooled_p = ks.sum() / ns.sum()
    observed_fracs = ks / ns
    observed_spread = observed_fracs.max() - observed_fracs.min()

    count_ge = 0
    for _ in range(n_boot):
        synth_fracs = np.array([monte_carlo_trials(pooled_p, int(n), rng).mean() for n in ns])
        synth_spread = synth_fracs.max() - synth_fracs.min()
        if synth_spread >= observed_spread:
            count_ge += 1
    p_value = count_ge / n_boot
    return pooled_p, observed_spread, p_value


def run_litcheck():
    """Sanity-check the Monte Carlo / probability layer against real (not synthetic)
    published repeat-trial nucleation counts. See the dataset constants above for
    citations and important caveats about what these numbers can and can't be used for."""
    print("=" * 78)
    print("LITERATURE REPEAT-TRIAL CHECK -- real published trial counts, not synthetic data.")
    print("These are NOT a fluence sweep and CANNOT calibrate kappa (see 'calibrate' mode")
    print("for that). This validates the Monte Carlo / Wilson-CI stochastic layer against")
    print("real reported nucleation statistics instead.")
    print("=" * 78)

    for label, dataset in [("KCl repeatability", KCL_REPEATABILITY_PMC10626568),
                            ("KBr additives", KBR_ADDITIVES_PMC10997761)]:
        print(f"\n--- {label} ---")
        print(f"Source: {dataset['citation']}")
        print(f"Note:   {dataset['note']}")
        for name, k, n in dataset["trials"]:
            lo, hi = wilson_ci(k, n)
            print(f"  {name:16s}  {k:3d}/{n:3d} = {k/n:.3f}   95% Wilson CI [{lo:.3f}, {hi:.3f}]")

    print("\n--- Binomial-consistency test (KCl repeatability only; pooling the KBr")
    print("    additive conditions would not be meaningful -- see note above) ---")
    pooled_p, observed_spread, p_value = bootstrap_spread_pvalue(
        KCL_REPEATABILITY_PMC10626568["trials"])
    print(f"Pooled probability across all three KCl sessions: {pooled_p:.3f}")
    print(f"Observed spread (max session fraction - min session fraction): {observed_spread:.3f}")
    print(f"Bootstrap p-value (chance pure binomial noise produces spread >= this): {p_value:.3f}")
    if p_value < 0.05:
        print("-> The real spread across sessions is larger than ordinary binomial sampling")
        print("   noise alone would typically produce. This points to genuine session-to-")
        print("   session physical variability (e.g. impurity content, alignment, ambient")
        print("   conditions) on top of the intrinsic stochasticity CNT already predicts --")
        print("   worth keeping in mind when you design your own repeat-trial protocol.")
    else:
        print("-> The real spread across sessions is consistent with ordinary binomial")
        print("   sampling noise around one shared probability -- no need to invoke anything")
        print("   beyond the stochastic nucleation model already implemented here.")


# =============================================================================
# Real S-dependent glycine NPLIN data (usable for an actual kappa fit)
# =============================================================================
# Javid, N.; Kendall, T.; Burns, I. S.; Sefcik, J. "Filtration Suppresses Laser-Induced
# Nucleation of Glycine in Aqueous Solutions." Cryst. Growth Des. 2016, 16 (8), 4196-4202.
# DOI: 10.1021/acs.cgd.6b00046 (values below extracted from the open-access author
# accepted manuscript, University of Strathclyde repository, Strathprints record 56853).
#
# CONFIDENCE LEVELS (read this before trusting any number below):
#   HIGH  -- stated explicitly, verbatim, in the paper's text/tables (laser wavelength,
#            power density, pulse duration, rep rate, exposure time, supersaturations,
#            Table 1 percentages, Table 2 fitted A/tau values and their 95% CIs).
#   LOW / ASSUMPTION -- NOT stated in the paper; estimated by us from vial geometry
#            described in the Methods section (11 mm vial diameter, 1 mm beam) to get an
#            illuminated volume, since the paper reports total vial volume (2 mL) but not
#            the illuminated sub-volume. Flagged explicitly wherever used.
#   PLACEHOLDER -- not measured in this paper or found elsewhere with confidence: solid-
#            liquid interfacial energy (gamma) and the CNT kinetic prefactor (A_prefactor)
#            for glycine. These remain the same kind of order-of-magnitude placeholder
#            used elsewhere in this script -- the fit below compensates for them by fitting
#            log(A) jointly with kappa, but that means the absolute fitted kappa is only as
#            trustworthy as that compensation, not an independently verified number.
#
# The paper models induction-time statistics as a mixture of two exponential nucleation
# regimes (Eq. 2): a fast, laser-associated regime with weight A and characteristic time
# tau1, and a slow, spontaneous regime with characteristic time tau2. The paper itself
# states the relation between characteristic time and nucleation rate: J*V = 1/tau (their
# Eq. following ref 27) -- i.e. exactly the same J*V*t Poisson-process picture this script
# already uses elsewhere. We treat "A" (the fitted fraction of samples attributable to the
# fast/laser regime) as the empirical laser-induced nucleation probability to fit against.
# CAVEAT: the paper's own preferred interpretation is that A reflects a heterogeneous
# sub-population of vials containing an activatable impurity, not a single homogeneous
# per-vial probability -- so treating A as "P(nucleate via light term)" for a homogeneous
# CNT+field model is an approximation, not a literal match to the authors' own mechanism.
GLYCINE_JAVID2016 = {
    "citation": ("Javid, Kendall, Burns & Sefcik, \"Filtration Suppresses Laser-Induced "
                 "Nucleation of Glycine in Aqueous Solutions,\" Cryst. Growth Des. 2016, "
                 "16, 4196-4202. DOI: 10.1021/acs.cgd.6b00046"),
    "laser": {
        "wavelength_nm": 1064.0,                 # HIGH confidence (stated)
        "power_density_W_m2": 0.47e9 * 1e4,       # HIGH confidence: 0.47 GW/cm^2 stated -> W/m^2
        "pulse_duration_s": 6e-9,                 # HIGH confidence (stated)
        "rep_rate_Hz": 10.0,                      # HIGH confidence (stated)
        "exposure_time_s": 60.0,                  # HIGH confidence: "1 minute of irradiation"
        "spot_diameter_m": 1e-3,                  # HIGH confidence: "1 mm circular pin-hole"
        "vial_diameter_m": 11e-3,                 # HIGH confidence (stated, for volume ASSUMPTION below)
    },
    "T_K": 298.15,  # HIGH confidence (stated: 298 K)
    # Table 1: total % crystallized within 5760 min (4 days) -- HIGH confidence.
    "table1_percent_crystallized": {
        1.4: {"non_filtered_irradiated": 21, "non_filtered_dark": 2},
        1.5: {"non_filtered_irradiated": 60, "non_filtered_dark": 11},
        1.6: {"non_filtered_irradiated": 94, "non_filtered_dark": 10},
    },
    # Table 2: biexponential fit to non-filtered irradiated samples -- HIGH confidence
    # (values and 95% CIs as printed in the paper). tau2 values marked (*) are flagged by
    # the AUTHORS as not statistically significant (insufficient data points) -- do not
    # use those two for anything quantitative.
    "table2_biexponential_fit": {
        1.4: {"A": 0.21, "A_ci95": 0.017, "td_min": 62.3, "tau1_min": 143.4, "tau1_ci95": 36.6,
              "tau2_min": 1.48e6, "tau2_significant": False},
        1.5: {"A": 0.45, "A_ci95": 0.012, "td_min": 37.4, "tau1_min": 117.0, "tau1_ci95": 8.6,
              "tau2_min": 1.14e4, "tau2_ci95": 0.16e4, "tau2_significant": True},
        1.6: {"A": 0.87, "A_ci95": 0.060, "td_min": 6.9, "tau1_min": 24.7, "tau1_ci95": 4.0,
              "tau2_min": 2.17e3, "tau2_significant": False},
    },
    # Table 3: S=1.5 only, spontaneous characteristic time tau_spont from non-irradiated
    # samples (single-exponential fit) -- HIGH confidence, only available at S=1.5.
    "tau_spont_min_S1p5": 3.43e4,
}


def estimate_illuminated_volume_m3(dataset=GLYCINE_JAVID2016):
    """ASSUMPTION, not stated in the paper: approximate illuminated volume as a cylinder of
    the stated beam diameter passing across the stated vial diameter. This is almost
    certainly an overestimate (real beam path/vial optics are more complex -- the paper
    itself notes the cylindrical vial acts as "a powerful cylindrical lens"), so treat this
    volume, and anything derived from it, as order-of-magnitude only."""
    d_beam = dataset["laser"]["spot_diameter_m"]
    path_length = dataset["laser"]["vial_diameter_m"]
    area = np.pi * (d_beam / 2) ** 2
    return area * path_length


# Sanity-check bound on nucleation rate ratio (fold-change). This is a heuristic, not a
# number taken from one specific paper: across the NPLIN literature cited in this script
# (Garetz/Myerson/Alexander/Camp/Knott and the reviews they wrote), reported laser-induced
# enhancements are consistently modest -- single-digit to low-hundreds fold changes in
# nucleation probability or rate, not many-orders-of-magnitude effects. 100x is deliberately
# at the upper end of that range (10-100x), not a loose/generous bound. A model result far
# outside this range is far more likely to be a fitting artifact than a real discovery.
# Adjust MAX_PLAUSIBLE_NPLIN_RATIO here if you have a better-justified bound for your system.
MAX_PLAUSIBLE_NPLIN_RATIO = 100.0

# The bound above must hold across this S range, not just near the real data points --
# otherwise a kappa that looks fine at S=1.4-1.7 could still blow up just outside it. S=1.0
# itself is a true mathematical singularity for ANY kappa>0 (zero supersaturation means zero
# thermodynamic driving force, so J0 -> 0 while the field term alone keeps J1 finite, making
# the ratio formally infinite) -- so we use S=1.05 as the practical near-threshold floor,
# not S=1.0 exactly. This is a genuine feature of the physics, not a way of dodging the check.
RATIO_CONSTRAINT_S_RANGE = (1.05, 2.0)


def rate_ratio_vs_kappa(S, kappa, material: MaterialParams, I_peak, n_solvent, n_solid):
    """J(with field)/J(no field) at a given S and kappa. Notably this does NOT depend on
    A_prefactor at all -- A is a common multiplicative factor in both J0 and J1 and cancels
    exactly in the ratio. Only gamma, v_m, T (material), S, and the field term (kappa,
    I_peak, n_solvent, n_solid) matter here."""
    u_field = field_energy_density(I_peak, n_solvent, n_solid, kappa)
    dG0 = delta_g_star(S, material.T, material.gamma, material.v_m)
    dG1 = delta_g_star_field(S, material.T, material.gamma, material.v_m, u_field)
    with np.errstate(over="ignore"):
        return np.exp(-(dG1 - dG0) / (K_B * material.T))


def find_kappa_for_ratio_cap(max_ratio, S_grid, material: MaterialParams, I_peak, n_solvent, n_solid):
    """Bisect for the largest kappa such that the ratio stays <= max_ratio across S_grid.
    The ratio is monotonically increasing in kappa (for n_solid > n_solvent), so this is a
    well-posed root-find: f(kappa) = max_over_S_grid(ratio) - max_ratio crosses zero once."""
    def f(kappa):
        ratios = [rate_ratio_vs_kappa(S, kappa, material, I_peak, n_solvent, n_solid) for S in S_grid]
        return max(ratios) - max_ratio

    kappa_hi = 1.0
    while f(kappa_hi) < 0 and kappa_hi < 1e12:
        kappa_hi *= 10
    if f(kappa_hi) < 0:
        return kappa_hi  # could not bracket a root; ratio never gets this large
    from scipy.optimize import brentq
    return brentq(f, 0.0, kappa_hi, xtol=1e-6 * kappa_hi if kappa_hi > 0 else 1e-6)


def run_calibrate_glycine():
    """Calibrate kappa (and the CNT prefactor A, jointly) against REAL published glycine
    NPLIN data (Javid et al. 2016 -- see GLYCINE_JAVID2016 above for full citation and
    confidence flags), replacing the earlier synthetic demo calibration. Then re-plots the
    nucleation-rate ratio using the fitted kappa and checks whether the predicted effect
    persists across the S=1.4-1.7 range where real glycine NPLIN experiments are conducted.

    IMPORTANT: an unconstrained fit of this model to only 3 real data points routinely runs
    away to an implausibly large kappa (see the explanation printed below) -- this function
    fits once unconstrained purely to show and diagnose that, then refits with an explicit
    plausibility bound on the resulting rate ratio and plots ONLY the bounded result."""
    ds = GLYCINE_JAVID2016
    print("=" * 78)
    print("CALIBRATING AGAINST REAL PUBLISHED DATA (glycine NPLIN)")
    print(f"Source: {ds['citation']}")
    print("This model is fit to real published NPLIN data on ORDINARY MOLECULAR CRYSTALS")
    print("(glycine) -- it is NOT fit or validated for diamond, carbon, or CO2 chemistry.")
    print("It tests plausibility of the general 'light biases nucleation' mechanism only.")
    print("=" * 78)

    laser = ds["laser"]
    fluence_Jcm2 = laser["power_density_W_m2"] * laser["pulse_duration_s"] / JCM2_TO_JM2

    print("\n--- Calibration data points used (S, fluence, nucleation probability) ---")
    print(f"{'S':>5} {'fluence (J/cm^2, per pulse)':>28} {'P (laser-induced fraction)':>28} {'95% CI':>10}  {'confidence':>10}")
    S_data, A_data, A_sigma = [], [], []
    for S, row in sorted(ds["table2_biexponential_fit"].items()):
        print(f"{S:5.2f} {fluence_Jcm2:28.4g} {row['A']:28.3f} {'+/-'+format(row['A_ci95'],'.3f'):>10}  {'REAL (Table 2)':>10}")
        S_data.append(S)
        A_data.append(row["A"])
        A_sigma.append(row["A_ci95"] / 1.96)  # 95% CI half-width -> 1-sigma
    S_data = np.array(S_data); A_data = np.array(A_data); A_sigma = np.array(A_sigma)
    print("All three (S, fluence, P) values above are REAL, stated explicitly in Javid et al.")
    print("2016 (Table 2 for P and its CI; Methods section for fluence, constant across all")
    print("three since the paper used a single fixed power density, not a fluence sweep).")
    print("Nothing in this table is digitized-from-a-figure or estimated -- the ONLY assumed")
    print("(not stated) quantity in this whole calibration is the illuminated volume below.")

    print("\n--- Laser parameters used (HIGH confidence, all stated explicitly in the paper) ---")
    print(f"wavelength={laser['wavelength_nm']:g} nm, power density={laser['power_density_W_m2']:.3g} W/m^2 "
          f"(0.47 GW/cm^2), pulse duration={laser['pulse_duration_s']:.3g} s, "
          f"rep rate={laser['rep_rate_Hz']:g} Hz, exposure={laser['exposure_time_s']:g} s "
          f"({int(laser['exposure_time_s']*laser['rep_rate_Hz'])} pulses)")

    V_assumed = estimate_illuminated_volume_m3(ds)
    print(f"\nILLUMINATED VOLUME: {V_assumed:.3e} m^3 ({V_assumed/UL_TO_M3:.3g} uL) "
          f"-- ASSUMPTION (not stated in paper), estimated from stated 1 mm beam diameter "
          f"across stated 11 mm vial diameter. Order-of-magnitude only.")

    # Material parameters: molar mass and a literature density value are real; interfacial
    # energy and CNT prefactor are still PLACEHOLDERS (see notes above) -- log(A) is fit
    # jointly with kappa below to compensate, so only the SHAPE of the fit (how well kappa
    # captures the S-dependence) should be trusted, not the absolute kappa magnitude.
    M_glycine = 75.07e-3  # kg/mol -- real, textbook value
    rho_glycine = 1.6e3   # kg/m^3 -- real, approximate literature value for glycine crystal density
    v_m = M_glycine / (rho_glycine * sc.N_A)
    gamma_placeholder = 1.3e-2  # J/m^2 -- PLACEHOLDER, not measured in this paper
    n_solvent = 1.33   # water, approximate
    n_solid = 1.55     # PLACEHOLDER -- typical organic-crystal order of magnitude, glycine value not confirmed

    material = MaterialParams(name="glycine (aqueous, Javid et al. 2016 conditions)",
                               gamma=gamma_placeholder, v_m=v_m, T=ds["T_K"],
                               A_prefactor=1e20, n_solvent=n_solvent, n_solid=n_solid)
    light = LightFieldParams(wavelength_nm=laser["wavelength_nm"],
                              pulse_energy_J=laser["power_density_W_m2"] * np.pi * (laser["spot_diameter_m"] / 2) ** 2 * laser["pulse_duration_s"],
                              spot_diameter_m=laser["spot_diameter_m"],
                              pulse_duration_s=laser["pulse_duration_s"],
                              n_pulses=int(laser["exposure_time_s"] * laser["rep_rate_Hz"]),
                              coupling_efficiency=1.0, pulse_shape_factor=1.0)
    geom = SampleGeometry(volume_m3=V_assumed, exposure_time_s=laser["exposure_time_s"])

    I_peak = laser["power_density_W_m2"]  # given directly, no need to derive from pulse energy

    def model(S_arr, kappa, log_A):
        A_pref = 10 ** log_A
        preds = []
        for S in S_arr:
            u_field = field_energy_density(I_peak, n_solvent, n_solid, kappa)
            J = nucleation_rate_field(S, material.T, material.gamma, material.v_m, A_pref, u_field)
            preds.append(probability_multi_pulse(J, geom.volume_m3, light.pulse_duration_s, light.n_pulses))
        return np.array(preds)

    # --- WHY kappa CAN RUN AWAY: kappa sits inside an exponential, twice removed ---
    # kappa enters linearly in u_field, which is ADDED linearly inside ln(S_eff). But
    # ln(S_eff) is then SQUARED in the denominator of dG* = 16*pi*gamma^3*v_m^2 / (3*(kT)^2*
    # ln(S_eff)^2), and dG* then sits inside J = A*exp(-dG*/kT), a NEGATIVE EXPONENTIAL. So
    # the chain is: kappa -(linear)-> ln(S_eff) -(inverse square)-> dG* -(exp)-> J. Near the
    # nucleation threshold (small ln S), that inverse-square step is extremely steep, so a
    # modest change in kappa can swing J (and hence the fitted probability) by many orders
    # of magnitude. With only 3 real data points and 2 free parameters (kappa, log10(A))
    # that are BOTH free to move along this steep, correlated surface -- and no independent
    # constraint on gamma or A for glycine -- curve_fit has more freedom than the data can
    # pin down. Left unconstrained (bounds up to +inf), it is free to run to whatever kappa
    # (paired with whatever A) minimizes the residual, even if that kappa implies an absurd
    # rate ratio. This is an ill-conditioned/underdetermined fit artifact, not a measurement.
    print("\n--- Why kappa can run away: it sits inside an exponential, twice removed ---")
    print("kappa -> (linear) -> ln(S_eff) -> (INVERSE SQUARE) -> dG* -> (NEGATIVE EXPONENTIAL) -> J")
    print("With only 3 data points and 2 free parameters (kappa, log10(A)) moving together on")
    print("this steep, correlated surface, and no independent constraint on gamma/A for glycine,")
    print("an unconstrained fit can drive kappa to a large value that reduces the residual without")
    print("that value being a meaningful physical measurement. Fitting this unconstrained first,")
    print("to show exactly that:")

    p0 = [1.0, 20.0]
    try:
        popt, pcov = curve_fit(model, S_data, A_data, p0=p0, sigma=A_sigma, absolute_sigma=True,
                                bounds=([0.0, 0.0], [np.inf, 40.0]), maxfev=20000)
        kappa_unc, logA_unc = popt
        kappa_unc_err = np.sqrt(pcov[0, 0])
        fit_ok = True
    except RuntimeError as e:
        print(f"\nFIT FAILED to converge: {e}")
        kappa_unc, logA_unc = p0
        kappa_unc_err = float("nan")
        fit_ok = False

    S_check_grid = np.linspace(*RATIO_CONSTRAINT_S_RANGE, 60)
    ratios_unc = [rate_ratio_vs_kappa(S, kappa_unc, material, I_peak, n_solvent, n_solid) for S in S_check_grid]
    max_ratio_unc = max(ratios_unc)
    print(f"\nUNCONSTRAINED fit: kappa = {kappa_unc:.4g} +/- {kappa_unc_err:.4g}")
    print(f"  -> implies a rate ratio up to {max_ratio_unc:.3g}x across S={RATIO_CONSTRAINT_S_RANGE[0]:g}-{RATIO_CONSTRAINT_S_RANGE[1]:g}")
    print(f"  -> sanity bound: real NPLIN literature enhancements are ~10-{MAX_PLAUSIBLE_NPLIN_RATIO:g}x, not orders")
    print(f"     of magnitude beyond that (see MAX_PLAUSIBLE_NPLIN_RATIO comment in the source).")
    print(f"  -> REJECTED as a fitting artifact (see explanation above), NOT plotted or used as a result.")

    # THE FIX: kappa_max is computed BEFORE fitting and passed directly into curve_fit's
    # `bounds=` argument, so the optimizer's search space itself never contains a kappa that
    # could produce a ratio above the plausibility bound anywhere in RATIO_CONSTRAINT_S_RANGE.
    # This is different from clipping the output after an unconstrained fit: if the true
    # least-squares optimum within [0, kappa_max] is AT the boundary kappa_max, that is a
    # real result of the constrained optimization (it means the data pulls harder than the
    # plausibility bound allows), not an artifact of clipping a larger unconstrained value.
    kappa_max = find_kappa_for_ratio_cap(MAX_PLAUSIBLE_NPLIN_RATIO, S_check_grid, material,
                                          I_peak, n_solvent, n_solid)
    print(f"\nConstrained search space: kappa in [0, {kappa_max:.4g}] -- the largest kappa for which")
    print(f"NO point in S={RATIO_CONSTRAINT_S_RANGE[0]:g}-{RATIO_CONSTRAINT_S_RANGE[1]:g} exceeds a "
          f"{MAX_PLAUSIBLE_NPLIN_RATIO:g}x ratio. Refitting kappa AND log10(A) jointly")
    print("with the optimizer restricted to this bound from the start (not fit-then-clip):")

    popt2, pcov2 = curve_fit(model, S_data, A_data, p0=[kappa_max * 0.5, 20.0], sigma=A_sigma,
                              absolute_sigma=True, bounds=([0.0, 0.0], [kappa_max, 40.0]), maxfev=20000)
    kappa_fit, logA_fit = popt2
    kappa_err, logA_err = np.sqrt(np.diag(pcov2))
    at_boundary = np.isclose(kappa_fit, kappa_max, rtol=1e-3)

    preds = model(S_data, kappa_fit, logA_fit)
    residuals = A_data - preds
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((A_data - A_data.mean())**2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    print("\n--- Constrained fit results (kappa search-bounded during optimization, not clipped after) ---")
    print(f"kappa      = {kappa_fit:.4g} +/- {kappa_err:.4g}" +
          ("  <- AT the constraint boundary" if at_boundary else "  (interior optimum, below the boundary)"))
    print(f"log10(A)   = {logA_fit:.4g} +/- {logA_err:.4g}  (A = {10**logA_fit:.3e} m^-3 s^-1, still a fitted placeholder)")
    print(f"\n{'S':>5} {'observed A':>12} {'model P':>12} {'residual':>10}")
    for S, obs, pred, res in zip(S_data, A_data, preds, residuals):
        print(f"{S:5.2f} {obs:12.3f} {pred:12.3f} {res:10.3f}")
    print(f"\nR^2 = {r_squared:.4f}  (only 3 data points / 2 free parameters -- 1 residual degree")
    print("of freedom, so treat this R^2 as descriptive, not a strong statistical validation.)")
    print(f"\nR^2 = {r_squared:.4f} vs. the REJECTED unconstrained fit's R^2 = "
          f"{1 - np.sum((A_data - model(S_data, kappa_unc, logA_unc))**2) / ss_tot:.4f} -- ", end="")
    if at_boundary:
        print("worse, as expected: the")
        print("optimizer wants to push kappa higher than the plausibility bound allows, and is")
        print("prevented from doing so. The optimum sitting exactly at the boundary is itself the")
        print("finding: fitting this model to this real data, under a literature-plausible ratio")
        print("bound, is fundamentally in tension -- the data wants more 'light effect' than a")
        print("physically modest coupling can honestly provide.")
        if r_squared < 0.5:
            print(f"\nFINDING: R^2={r_squared:.2f} under the {MAX_PLAUSIBLE_NPLIN_RATIO:g}x-bounded constraint is a poor fit.")
            print("This says the simple kappa-only perturbative light-field term, held to a")
            print("literature-plausible rate-ratio bound, CANNOT reasonably reproduce the observed")
            print("S-dependence in this real glycine dataset (P rising from 0.21 to 0.87 across")
            print("S=1.4-1.6). That is itself a real result, not a failure of this exercise: it")
            print("suggests either (a) the true mechanism behind these numbers is not a small,")
            print("uniform free-energy perturbation of this form, or (b) the placeholder gamma/A")
            print("values are far enough from glycine's real values that no plausible kappa can")
            print("compensate. Recall the source paper's OWN conclusion (filtration suppresses the")
            print("effect) -- pointing at an impurity/thermocavitation-mediated mechanism rather")
            print("than direct field-cluster coupling -- which this script does not model at all.")
    else:
        magnitude_change = abs(np.log10(max(kappa_unc, 1e-12)) - np.log10(max(kappa_fit, 1e-12)))
        print("comparable -- the plausibility")
        print("bound was NOT actually binding; the interior optimum (kappa below the boundary)")
        print("already satisfied it.")
        if magnitude_change >= 1.0:
            print(f"kappa changed by about {magnitude_change:.1f} orders of magnitude between the")
            print(f"unconstrained ({kappa_unc:.4g}) and constrained ({kappa_fit:.4g}) fits, yet R^2")
            print("barely moved. That means log10(A) (the free CNT prefactor), NOT kappa, is doing")
            print("almost all the work of matching the observed S-dependence: ordinary CNT's own")
            print("(ln S)^-2 sensitivity, with little-to-no light-field contribution, already")
            print("reproduces most of the shape of this real data. See the perturbation check below")
            print("for how small the light term's actual contribution turns out to be.")
        else:
            print(f"kappa changed only modestly between the unconstrained ({kappa_unc:.4g}) and")
            print(f"constrained ({kappa_fit:.4g}) fits -- the plausibility bound was not very")
            print("restrictive for this functional form: the 'runaway' unconstrained value was")
            print("already close to a physically-plausible magnitude here, unlike cases where an")
            print("unconstrained fit runs away by many orders of magnitude. Check the perturbation")
            print("check below to judge whether this fitted kappa is a small or large perturbation")
            print("relative to the thermodynamic driving force.")

    # Ratio plot with the properly constrained (not clipped) kappa, over the full range the
    # constraint was checked against, with the real experimental window highlighted.
    light.coupling_efficiency = kappa_fit
    material.A_prefactor = 10 ** logA_fit
    disclaimer = (f"Fit to real glycine data only -- NOT diamond/carbon/CO2; kappa search-bounded "
                   f"to keep ratio <= {MAX_PLAUSIBLE_NPLIN_RATIO:g}x across S={RATIO_CONSTRAINT_S_RANGE[0]:g}-{RATIO_CONSTRAINT_S_RANGE[1]:g}")
    S_range, J0, J1 = report_and_plot_rate_vs_S(
        material, light, geom, "glycine_calibrated_ratio.png",
        S_min=RATIO_CONSTRAINT_S_RANGE[0], S_max=RATIO_CONSTRAINT_S_RANGE[1],
        highlight_S_range=(1.4, 1.7), highlight_label="target range S=1.4-1.7",
        title_disclaimer=disclaimer,
    )

    # CRITICAL PHYSICAL-PLAUSIBILITY CHECK, separate from the ratio magnitude itself: is the
    # fitted kappa perturbing the free energy a little (as the "light nudges nucleation"
    # framing assumes), or is it required to be as large as the thermodynamic driving force
    # itself? If the latter, the fit "succeeding" mostly reflects the placeholder gamma/A
    # values being compensated for, not a physically modest light-coupling effect.
    u_field_fit = field_energy_density(I_peak, n_solvent, n_solid, kappa_fit)
    rel_shift_fit = u_field_fit * material.v_m / (K_B * material.T)
    lnS_mid = np.log(1.5)
    shift_fraction = rel_shift_fit / lnS_mid
    print("\n--- Is the fitted kappa a small perturbation, or does it dominate? ---")
    print(f"Light-induced shift to ln(S): {rel_shift_fit:.4g}  vs.  ln(1.5) = {lnS_mid:.4g}")
    print(f"Ratio (shift / ln(S)): {shift_fraction:.3g}")
    if shift_fraction > 0.3:
        print("This is NOT a small perturbation -- the fitted light term is comparable to, or")
        print("larger than, the thermodynamic driving force itself. That means this fit result")
        print("is dominated by compensating for the placeholder gamma/A_prefactor values (which")
        print("were NOT independently measured for glycine in this exercise), not evidence of a")
        print("physically modest 'light nudges nucleation' effect. Treat the PASS/FAIL verdict")
        print("below with this firmly in mind.")
    else:
        print("The fitted light term is a modest perturbation relative to the thermodynamic")
        print("driving force -- consistent with the 'light nudges nucleation' framing this")
        print("model was built around.")

    print("\n" + "=" * 78)
    print("PASS/FAIL CHECK: does the fitted-kappa effect persist across S=1.4-1.7,")
    print("the range where real glycine NPLIN experiments are actually conducted?")
    print("=" * 78)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_full = J1 / J0
    mask_1417 = (S_range >= 1.4) & (S_range <= 1.7)
    if mask_1417.any():
        r_lo, r_hi = ratio_full[mask_1417].min(), ratio_full[mask_1417].max()
        dev_pct = max(abs(r_lo - 1), abs(r_hi - 1)) * 100
        if dev_pct >= 5.0:
            verdict = "PASS"
        elif dev_pct >= 1.0:
            verdict = "MARGINAL"
        else:
            verdict = "FAIL"
        print(f"Ratio within S=1.4-1.7: {r_lo:.4g} to {r_hi:.4g} (max deviation from 1: {dev_pct:.3g}%)")
        if verdict == "PASS" and shift_fraction > 0.3:
            print("Verdict: PASS (ratio far from 1) BUT NOT PHYSICALLY MODEST -- ", end="")
            print("see the perturbation check above: this large ratio is achieved by making the")
            print("light term comparable to the thermodynamic driving force itself, not by a")
            print("small nudge. Read this as 'the model *can* be forced to fit the data', not as")
            print("'a modest light-coupling hypothesis is confirmed'.")
            verdict = "PASS_NOT_MODEST"
        elif verdict == "PASS":
            print(f"Verdict: {verdict} -- ", end="")
            print("the fitted effect is meaningfully present across the real experimental range.")
        elif verdict == "MARGINAL":
            print(f"Verdict: {verdict} -- ", end="")
            print("a small (1-5%) effect persists -- probably too small to reliably distinguish")
            print("from trial-to-trial binomial noise without a large number of repeat trials.")
        else:
            print(f"Verdict: {verdict} -- ", end="")
            print("the fitted effect has decayed to <1% across the real experimental range --")
            print("consistent with the earlier finding that this mechanism's leverage is")
            print("concentrated in a narrow band near the classical nucleation threshold, not")
            print("across the S range real NPLIN experiments actually use.")
    else:
        verdict = "N/A"
        print("S=1.4-1.7 falls outside the plotted S range -- re-run with a wider S_max.")

    print("\n" + "=" * 78)
    print("PLAIN-LANGUAGE SUMMARY")
    print("=" * 78)
    print("This is a model fit to real published glycine (ordinary molecular crystal) NPLIN")
    print("data, with only 3 real data points and 2 free parameters -- a weak statistical")
    print("test, not a strong validation. Within that limitation:")
    if not fit_ok:
        print("- The fit did not converge; no verdict available (see FIT FAILED message above).")
    elif verdict == "PASS_NOT_MODEST":
        print("- HONEST ANSWER: this is inconclusive-to-negative, not a confirmation. The fit CAN")
        print("  reproduce the real S-dependence, but only by making the light term as large as")
        print("  the thermodynamic driving force itself -- which is not what 'light nudges")
        print("  nucleation' is supposed to mean, and mostly reflects the placeholder gamma/A")
        print("  values (not independently measured here) absorbing the mismatch. This does NOT")
        print("  support the hypothesis in the modest-perturbation sense it was proposed in.")
    elif verdict == "PASS":
        print("- The calibrated model DOES support the hypothesis that a light-field term can")
        print("  meaningfully bias nucleation rate at realistic glycine NPLIN supersaturations.")
    elif verdict == "MARGINAL":
        print("- The calibrated model gives WEAK, inconclusive support: a real but small effect")
        print("  at realistic supersaturations, likely hard to distinguish from noise experimentally.")
    elif verdict == "FAIL":
        print("- The calibrated model DOES NOT support a meaningful effect at realistic glycine")
        print("  NPLIN supersaturations -- fitting kappa to reproduce the observed S-dependence in")
        print("  Table 2 pulls the model's sensitive region towards the data, but that sensitive")
        print("  region is still concentrated near the nucleation threshold, not spread across")
        print("  S=1.4-1.7.")
    else:
        print("- S=1.4-1.7 was outside the plotted range; no verdict available.")
    print("This says nothing about diamond, carbon, or CO2 systems -- only about whether the")
    print("general mechanism, calibrated against real molecular-crystal data, is plausible in")
    print("its own domain.")

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_curve = J1 / J0
    return {"label": "Original (barrier-modifying, fixed coupling)", "kappa": kappa_fit,
            "logA": logA_fit, "r_squared": r_squared, "verdict": verdict,
            "S_range": S_range, "ratio": ratio_curve}


# =============================================================================
# Calibration against real data
# =============================================================================
def fit_calibration(fluence_arr_Jm2, prob_arr, S, material: MaterialParams,
                     light_template: LightFieldParams, geom: SampleGeometry,
                     fit_params=("kappa",), p0=None, bounds=None):
    """Fit free parameters (subset of {'kappa','logA'}) to observed nucleation
    probability vs. fluence data at a fixed, known supersaturation S. Parameters
    not being fit are held at the values in `material`/`light_template`.

    Returns (fitted_dict, pcov, model_fn).
    """
    param_names = list(fit_params)

    def unpack(theta):
        d = dict(zip(param_names, theta))
        kappa = d.get("kappa", light_template.coupling_efficiency)
        A = 10 ** d["logA"] if "logA" in d else material.A_prefactor
        return kappa, A

    def model(H, *theta):
        kappa, A = unpack(theta)
        I_peak = peak_intensity_from_fluence(H, light_template.pulse_duration_s,
                                              light_template.pulse_shape_factor)
        u_field = field_energy_density(I_peak, material.n_solvent, material.n_solid, kappa)
        J = nucleation_rate_field(S, material.T, material.gamma, material.v_m, A, u_field)
        return probability_multi_pulse(J, geom.volume_m3, light_template.pulse_duration_s,
                                        light_template.n_pulses)

    if p0 is None:
        p0 = [0.5 if p == "kappa" else np.log10(material.A_prefactor) for p in param_names]
    if bounds is None:
        lo = [0.0 if p == "kappa" else -np.inf for p in param_names]
        hi = [np.inf for _ in param_names]
        bounds = (lo, hi)

    popt, pcov = curve_fit(model, fluence_arr_Jm2, prob_arr, p0=p0, bounds=bounds, maxfev=20000)
    return dict(zip(param_names, popt)), pcov, model


def load_calibration_csv(path):
    """CSV columns: fluence_J_cm2, probability, [n_trials]. Returns fluence in J/m^2."""
    fluences, probs, n_trials = [], [], []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            fluences.append(float(row["fluence_J_cm2"]) * JCM2_TO_JM2)
            probs.append(float(row["probability"]))
            n_trials.append(float(row["n_trials"]) if "n_trials" in row and row["n_trials"] not in (None, "") else np.nan)
    return np.array(fluences), np.array(probs), np.array(n_trials)


# =============================================================================
# Interactive prompts
# =============================================================================
def ask_float(prompt_text, default=None, source=None):
    suffix = ""
    if default is not None:
        suffix = f" [default {default:g}"
        if source:
            suffix += f" -- {source}"
        suffix += "]"
    while True:
        raw = input(f"{prompt_text}{suffix}: ").strip()
        if raw == "" and default is not None:
            return default
        try:
            return float(raw)
        except ValueError:
            print("  Please enter a numeric value.")


def ask_int(prompt_text, default=None):
    suffix = f" [default {default}]" if default is not None else ""
    while True:
        raw = input(f"{prompt_text}{suffix}: ").strip()
        if raw == "" and default is not None:
            return default
        try:
            return int(raw)
        except ValueError:
            print("  Please enter an integer.")


def ask_str(prompt_text, default=None):
    suffix = f" [default '{default}']" if default is not None else ""
    raw = input(f"{prompt_text}{suffix}: ").strip()
    return raw if raw else default


def ask_yes_no(prompt_text, default=False):
    suffix = " [Y/n]" if default else " [y/N]"
    raw = input(f"{prompt_text}{suffix}: ").strip().lower()
    if raw == "":
        return default
    return raw.startswith("y")


def gather_params_interactively():
    print("\n=== Material / solution parameters ===")
    print("(These are compound- and polymorph-specific. Use your own measured")
    print(" or literature-sourced values -- defaults offered here are ONLY")
    print(" illustrative placeholders, clearly marked as such.)\n")

    name = ask_str("Compound / system name (e.g. 'glycine, aqueous, alpha polymorph')",
                    default="unnamed system")

    gamma = ask_float(
        "Solid-liquid interfacial free energy, gamma (J/m^2)",
        default=2e-3,
        source=("PLACEHOLDER ONLY: order of magnitude typical of small organic "
                "molecules crystallizing from aqueous solution (general texts, e.g. "
                "Mullin, Crystallization). Replace with a value reported/measured "
                "for YOUR compound, polymorph, and solvent."),
    )

    have_vm = ask_yes_no("Do you already know the molecular volume in the solid, v_m, directly?",
                          default=False)
    if have_vm:
        v_m = ask_float("Molecular volume in the solid, v_m (m^3/molecule)", default=8e-29,
                         source="PLACEHOLDER -- typical order of magnitude for a small organic molecule crystal.")
    else:
        M = ask_float("Molar mass of the crystallizing species (g/mol)", default=75.07,
                       source="glycine molar mass, for reference only -- use your compound's value")
        rho = ask_float("Crystal density (g/cm^3)", default=1.43,
                         source="PLACEHOLDER order of magnitude -- look up your specific polymorph's density")
        v_m = (M * 1e-3) / (rho * 1e3) / sc.N_A  # m^3/molecule
        print(f"  -> computed v_m = {v_m:.3e} m^3/molecule")

    T = ask_float("Temperature, T (K)", default=298.15, source="room temperature")

    A_prefactor = ask_float(
        "CNT kinetic prefactor, A (1/(m^3 s))",
        default=1e30,
        source=("PLACEHOLDER: textbook order-of-magnitude range is roughly 1e25-1e35 "
                "m^-3 s^-1 (Kashchiev, Nucleation, 2000; Mullin, Crystallization, 2001). "
                "This is the single most uncertain CNT parameter -- strongly prefer "
                "fitting it (via 'calibrate' mode) to your own induction-time or "
                "nucleation-probability data rather than trusting this default."),
    )

    n_solvent = ask_float("Refractive index of the solvent/mother liquor at the laser wavelength",
                           default=1.33, source="approx. water in the visible; use your solvent's measured value")
    n_solid = ask_float("Refractive index of the crystal (cluster) at the laser wavelength",
                         default=1.53, source="PLACEHOLDER order of magnitude for a small organic molecular crystal")

    print("\n=== Light-field (laser) parameters ===\n")
    wavelength = ask_float("Laser wavelength (nm)", default=532.0)
    pulse_energy_mJ = ask_float("Pulse energy (mJ)", default=10.0)
    spot_diam_mm = ask_float("Beam spot diameter at the sample (mm)", default=1.0)
    pulse_duration_ns = ask_float("Pulse duration (ns)", default=6.0,
                                   source="typical Q-switched Nd:YAG pulse width")
    n_pulses = ask_int("Number of pulses applied per trial", default=1)
    kappa = ask_float(
        "Light-coupling efficiency, kappa (dimensionless, tunable)",
        default=0.1,
        source=("This is YOUR hypothesis parameter, not a literature constant: 0 = no "
                "effect (ordinary CNT), 1 = take the nominal dielectric-sphere field "
                "energy at face value. The underlying mechanism (isotropic electronic "
                "polarization / 'DP model') is proposed in Alexander & Camp, J. Chem. "
                "Phys. 150, 040901 (2019). There is no settled consensus value -- use "
                "'calibrate' mode to fit kappa against real data."),
    )
    shape_factor = ask_float("Temporal pulse-shape peak-intensity correction factor", default=1.0,
                              source="1.0 = flat-top approx; ~1.1-1.2 if you want to approximate a Gaussian-in-time pulse's higher peak vs. its average")

    print("\n=== Sample / observation geometry ===\n")
    volume_uL = ask_float("Illuminated / observed solution volume (microliters)", default=100.0)
    exposure_s = ask_float("Observation/holding time per trial after the pulse(s) (s)", default=1.0)

    material = MaterialParams(name=name, gamma=gamma, v_m=v_m, T=T, A_prefactor=A_prefactor,
                               n_solvent=n_solvent, n_solid=n_solid)
    light = LightFieldParams(wavelength_nm=wavelength, pulse_energy_J=pulse_energy_mJ * MJ_TO_J,
                              spot_diameter_m=spot_diam_mm * MM_TO_M,
                              pulse_duration_s=pulse_duration_ns * NS_TO_S, n_pulses=n_pulses,
                              coupling_efficiency=kappa, pulse_shape_factor=shape_factor)
    geom = SampleGeometry(volume_m3=volume_uL * UL_TO_M3, exposure_time_s=exposure_s)
    return material, light, geom


def load_config(path):
    data = json.loads(Path(path).read_text())
    material = MaterialParams(**data["material"])
    light = LightFieldParams(**data["light"])
    geom = SampleGeometry(**data["geometry"])
    return material, light, geom


def save_config(path, material, light, geom):
    data = {"material": asdict(material), "light": asdict(light), "geometry": asdict(geom)}
    Path(path).write_text(json.dumps(data, indent=2))
    print(f"Saved parameters to {path}")


def get_params(args):
    if getattr(args, "config", None) and Path(args.config).exists():
        print(f"Loading parameters from {args.config}")
        material, light, geom = load_config(args.config)
    else:
        material, light, geom = gather_params_interactively()
    if getattr(args, "save_config", None):
        save_config(args.save_config, material, light, geom)
    return material, light, geom


# =============================================================================
# Plotting / reporting
# =============================================================================
def report_and_plot_rate_vs_S(material: MaterialParams, light: LightFieldParams,
                               geom: SampleGeometry, out_path: str, S_max=4.0, S_min=1.01, n_points=300,
                               highlight_S_range=None, highlight_label=None, title_disclaimer=None):
    S_range = np.linspace(S_min, S_max, n_points)

    J0 = nucleation_rate(S_range, material.T, material.gamma, material.v_m, material.A_prefactor)
    u_field = u_field_from_params(material, light)
    J1 = nucleation_rate_field(S_range, material.T, material.gamma, material.v_m,
                                material.A_prefactor, u_field)

    H = fluence_from_pulse(light.pulse_energy_J, light.spot_diameter_m)
    I_peak = peak_intensity_from_fluence(H, light.pulse_duration_s, light.pulse_shape_factor)

    # Ratio (fold-change) plot: this is the number that answers "does kappa do anything
    # meaningful", since J0 and J1 plotted separately on a log scale are often visually
    # indistinguishable even when nowhere near equal, and can also look "different" when
    # the actual ratio is 1.0000006 -- i.e. both directions of misleading impression are
    # possible from the raw-curve overlay. The ratio strips that ambiguity out directly.
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = J1 / J0
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(S_range, ratio, lw=2, color="tab:purple",
                label=f"J(with field) / J(no field), kappa={light.coupling_efficiency:g}")
    ax.axhline(1.0, color="k", ls="--", lw=1, label="no effect (ratio = 1)")
    if highlight_S_range is not None:
        ax.axvspan(highlight_S_range[0], highlight_S_range[1], color="tab:green", alpha=0.15,
                   label=highlight_label or f"real experimental range S={highlight_S_range[0]:g}-{highlight_S_range[1]:g}")
    ax.set_xlabel("Supersaturation, S")
    ax.set_ylabel("Nucleation rate ratio, J(with field) / J(no field)")
    title = f"{material.name}\nfluence={H/JCM2_TO_JM2:.3g} J/cm^2, I_peak={I_peak:.3g} W/m^2, wavelength={light.wavelength_nm:g} nm"
    if title_disclaimer:
        title += f"\n{title_disclaimer}"
    ax.set_title(title, fontsize=9, wrap=True)
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    max_ratio, min_ratio = np.nanmax(ratio), np.nanmin(ratio)
    S_at_max = S_range[np.nanargmax(ratio)]
    max_dev_pct = max(abs(max_ratio - 1), abs(min_ratio - 1)) * 100
    ax.text(0.02, 0.02,
            f"ratio range over this S sweep: {min_ratio:.6g} - {max_ratio:.6g}\n"
            f"peak ratio {max_ratio:.6g} occurs at S={S_at_max:.4g}\n"
            f"(largest deviation from 1: {max_dev_pct:.4g}%)",
            transform=ax.transAxes, fontsize=8, va="bottom",
            bbox=dict(boxstyle="round", fc="white", alpha=0.8))

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"\nSaved plot to {out_path}")
    print(f"Ratio range over this S sweep: {min_ratio:.6g} to {max_ratio:.6g} "
          f"(largest deviation from 1: {max_dev_pct:.4g}%)")
    print(f"Peak ratio {max_ratio:.6g} occurs at S={S_at_max:.4g}")
    if highlight_S_range is not None:
        in_range_mask = (S_range >= highlight_S_range[0]) & (S_range <= highlight_S_range[1])
        if in_range_mask.any():
            range_max_dev_pct = max(abs(ratio[in_range_mask].max() - 1), abs(ratio[in_range_mask].min() - 1)) * 100
            print(f"Within the real experimental range S={highlight_S_range[0]:g}-{highlight_S_range[1]:g}: "
                  f"ratio spans {ratio[in_range_mask].min():.6g} to {ratio[in_range_mask].max():.6g} "
                  f"(largest deviation from 1 in-range: {range_max_dev_pct:.4g}%)")
    if max_dev_pct < 1.0:
        print("This is a sub-1% effect across the whole sweep at this kappa/fluence -- i.e.")
        print("kappa would need to be substantially larger (or the fluence substantially")
        print("higher) to produce anything you'd expect to resolve against trial-to-trial")
        print("nucleation-probability noise in a real experiment.")
    try:
        plt.show()
    except Exception:
        pass

    print("\n=== Fold-change in nucleation rate at selected supersaturations ===")
    print(f"{'S':>8} {'J (no field)':>16} {'J (with field)':>16} {'fold-change':>14} "
          f"{'r* no field (nm)':>18} {'r* with field (nm)':>20}")
    for S_check in [1.2, 1.5, 2.0, 3.0]:
        if S_check > S_max:
            continue
        j0 = nucleation_rate(S_check, material.T, material.gamma, material.v_m, material.A_prefactor)
        j1 = nucleation_rate_field(S_check, material.T, material.gamma, material.v_m,
                                    material.A_prefactor, u_field)
        r0 = critical_radius(S_check, material.T, material.gamma, material.v_m) * 1e9
        r1 = critical_radius_field(S_check, material.T, material.gamma, material.v_m, u_field) * 1e9
        fold = j1 / j0 if j0 > 0 else np.nan
        print(f"{S_check:8.2f} {j0:16.3e} {j1:16.3e} {fold:14.3e} {r0:18.3f} {r1:20.3f}")

    rel_shift = abs(u_field) * material.v_m / (K_B * material.T)
    print(f"\nu_field = {u_field:.3e} J/m^3 "
          f"({'promotes' if u_field > 0 else 'suppresses'} nucleation at this fluence)")
    print(f"This shifts ln(S) by {rel_shift:.3e} (kappa={light.coupling_efficiency:g}). "
          f"For reference, ln(S) itself is order 0.01-1 across typical supersaturations --")
    print(f"if rel_shift is many orders of magnitude smaller than ln(S), the nominal effect")
    print(f"at these pulse parameters is negligible regardless of the fold-change numbers above")
    print(f"looking non-trivial (rounding can hide a genuinely tiny effect). Compare rel_shift")
    print(f"to your actual ln(S) working point before concluding the light field matters here.")
    return S_range, J0, J1


def run_demo():
    print(__doc__.split("USAGE")[0])
    print("=" * 78)
    print("DEMO MODE: everything below uses CLEARLY LABELED placeholder parameters")
    print("and SYNTHETIC (computer-generated, not real) calibration data. This is")
    print("only meant to exercise the full pipeline end to end.")
    print("=" * 78)

    material = MaterialParams(
        name="illustrative glycine-like system (DEMO / not validated)",
        gamma=1.3e-2, v_m=8.0e-29, T=298.15, A_prefactor=1e18,
        n_solvent=1.33, n_solid=1.53,
    )
    light = LightFieldParams(
        wavelength_nm=532.0, pulse_energy_J=10e-3, spot_diameter_m=1e-3,
        pulse_duration_s=6e-9, n_pulses=1, coupling_efficiency=0.15, pulse_shape_factor=1.0,
    )
    geom = SampleGeometry(volume_m3=100e-9, exposure_time_s=1.0)

    report_and_plot_rate_vs_S(material, light, geom, "demo_rate_vs_S.png")

    print("\n" + "=" * 78)
    print("NOTE on the fold-change table just printed: at the realistic pulse energy")
    print("used above (10 mJ, 1 mm spot, 6 ns, kappa=0.15), the nominal DP-model field")
    print("perturbation is many orders of magnitude smaller than ln(S) itself (see the")
    print("rel_shift line above) -- i.e. at achievable unfocused-beam pulse energies this")
    print("mechanism, taken at face value, predicts an experimentally invisible effect.")
    print("This mirrors the real literature: Knott et al., J. Chem. Phys. 134, 154501")
    print("(2011) found the (related) optical-Kerr mechanism needs field strengths far")
    print("above experimental values to matter. Getting a VISIBLE effect out of this")
    print("model requires either much higher intensity (tighter focus / more energy --")
    print("watch out for optical breakdown/cavitation once you do that, which is a")
    print("different, unmodeled mechanism) or kappa >> 1, i.e. assuming the real")
    print("coupling is far stronger than the naive dielectric-sphere estimate.")
    print("=" * 78)

    print("\n" + "=" * 78)
    print("Now running the REAL calibration (Javid et al. 2016 glycine data) instead of a")
    print("synthetic self-test. This is the exact same calibration 'calibrate-glycine' mode")
    print("runs on its own -- demo mode just chains it here so the whole pipeline (including")
    print("the fitting machinery) is exercised end to end against real data by default.")
    print("=" * 78)

    run_calibrate_glycine()

    print("\nTo calibrate against your OWN data instead: 'python nplin_cnt_model.py calibrate")
    print("--data yourfile.csv' with a CSV of columns: fluence_J_cm2,probability[,n_trials]")


def run_explore(args):
    material, light, geom = get_params(args)
    report_and_plot_rate_vs_S(material, light, geom, args.out, S_max=args.s_max)


def run_calibrate(args):
    material, light, geom = get_params(args)
    S = ask_float("Supersaturation S at which this calibration dataset was measured", default=None)
    fluence_arr, prob_arr, n_trials_arr = load_calibration_csv(args.data)

    fit_params = tuple(p.strip() for p in args.fit.split(","))
    fitted, pcov, model_fn = fit_calibration(fluence_arr, prob_arr, S, material, light, geom,
                                              fit_params=fit_params)

    print("\n=== Calibration fit results ===")
    for i, name in enumerate(fit_params):
        err = np.sqrt(pcov[i, i])
        print(f"  {name} = {fitted[name]:.4g} +/- {err:.4g}")

    if "kappa" in fitted:
        light.coupling_efficiency = fitted["kappa"]
    if "logA" in fitted:
        material.A_prefactor = 10 ** fitted["logA"]
    if getattr(args, "save_config", None):
        save_config(args.save_config, material, light, geom)

    Hfine = np.linspace(fluence_arr.min(), fluence_arr.max(), 200)
    fig, ax = plt.subplots(figsize=(7, 5))
    if np.all(~np.isnan(n_trials_arr)):
        yerr_lo, yerr_hi = [], []
        for p, n in zip(prob_arr, n_trials_arr):
            lo, hi = wilson_ci(round(p * n), n)
            yerr_lo.append(p - lo); yerr_hi.append(hi - p)
        ax.errorbar(fluence_arr / JCM2_TO_JM2, prob_arr, yerr=[yerr_lo, yerr_hi], fmt="o", label="Data")
    else:
        ax.plot(fluence_arr / JCM2_TO_JM2, prob_arr, "o", label="Data")
    ax.plot(Hfine / JCM2_TO_JM2, model_fn(Hfine, *[fitted[p] for p in fit_params]), "-", label="Fit")
    ax.set_xlabel("Laser fluence (J/cm^2)")
    ax.set_ylabel("Nucleation probability")
    ax.set_title(f"Calibration fit: {material.name} at S={S:g}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Saved plot to {args.out}")
    try:
        plt.show()
    except Exception:
        pass


def run_montecarlo(args):
    material, light, geom = get_params(args)
    S_range = np.linspace(1.05, args.s_max, 12)
    u_field = u_field_from_params(material, light)

    rng = np.random.default_rng()
    P0_list, P1_list = [], []
    obs0, obs1 = [], []
    ci0_lo, ci0_hi, ci1_lo, ci1_hi = [], [], [], []

    for S in S_range:
        J0 = nucleation_rate(S, material.T, material.gamma, material.v_m, material.A_prefactor)
        J1 = nucleation_rate_field(S, material.T, material.gamma, material.v_m,
                                    material.A_prefactor, u_field)
        P0 = probability_multi_pulse(J0, geom.volume_m3, light.pulse_duration_s, light.n_pulses) \
            if light.n_pulses > 0 else nucleation_probability(J0, geom.volume_m3, geom.exposure_time_s)
        P1 = probability_multi_pulse(J1, geom.volume_m3, light.pulse_duration_s, light.n_pulses) \
            if light.n_pulses > 0 else nucleation_probability(J1, geom.volume_m3, geom.exposure_time_s)
        P0_list.append(P0); P1_list.append(P1)

        t0 = monte_carlo_trials(P0, args.n_trials, rng)
        t1 = monte_carlo_trials(P1, args.n_trials, rng)
        k0, k1 = t0.sum(), t1.sum()
        obs0.append(k0 / args.n_trials); obs1.append(k1 / args.n_trials)
        lo, hi = wilson_ci(k0, args.n_trials); ci0_lo.append(lo); ci0_hi.append(hi)
        lo, hi = wilson_ci(k1, args.n_trials); ci1_lo.append(lo); ci1_hi.append(hi)

    P0_list = np.array(P0_list); P1_list = np.array(P1_list)
    obs0 = np.array(obs0); obs1 = np.array(obs1)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(S_range, P0_list, "b-", label="Analytic P, no field")
    ax.plot(S_range, P1_list, "r-", label="Analytic P, with field")
    ax.errorbar(S_range, obs0, yerr=[obs0 - np.array(ci0_lo), np.array(ci0_hi) - obs0],
                fmt="bo", alpha=0.6, label=f"MC observed fraction, no field (n={args.n_trials}/pt)")
    ax.errorbar(S_range, obs1, yerr=[obs1 - np.array(ci1_lo), np.array(ci1_hi) - obs1],
                fmt="ro", alpha=0.6, label=f"MC observed fraction, with field (n={args.n_trials}/pt)")
    ax.set_xlabel("Supersaturation, S")
    ax.set_ylabel("Nucleation probability")
    ax.set_title(f"{material.name}: stochastic trials vs. analytic probability")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Saved plot to {args.out}")
    try:
        plt.show()
    except Exception:
        pass


# =============================================================================
# CLI
# =============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="CNT + light-field-term NPLIN plausibility simulator (see module docstring for physics/refs).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="mode", required=True)

    sub.add_parser("demo", help="Self-contained illustrative demo with placeholder params and synthetic data.")

    p_exp = sub.add_parser("explore", help="Enter your system's parameters; plot J vs S with/without the light field.")
    p_exp.add_argument("--config", type=str, default=None)
    p_exp.add_argument("--save-config", type=str, default=None)
    p_exp.add_argument("--out", type=str, default="nplin_rate_vs_S.png")
    p_exp.add_argument("--s-max", type=float, default=4.0)

    p_cal = sub.add_parser("calibrate", help="Fit the light-field coupling to real nucleation-probability-vs-fluence data.")
    p_cal.add_argument("--data", type=str, required=True,
                        help="CSV with header: fluence_J_cm2,probability[,n_trials]")
    p_cal.add_argument("--config", type=str, default=None)
    p_cal.add_argument("--save-config", type=str, default=None)
    p_cal.add_argument("--fit", type=str, default="kappa",
                        help="Comma-separated subset of {kappa,logA} to fit (default: kappa)")
    p_cal.add_argument("--out", type=str, default="nplin_calibration_fit.png")

    p_mc = sub.add_parser("montecarlo", help="Stochastic trials vs. analytic probability, swept over supersaturation.")
    p_mc.add_argument("--config", type=str, default=None)
    p_mc.add_argument("--save-config", type=str, default=None)
    p_mc.add_argument("--n-trials", type=int, default=50)
    p_mc.add_argument("--s-max", type=float, default=4.0)
    p_mc.add_argument("--out", type=str, default="nplin_montecarlo.png")

    sub.add_parser("litcheck", help="Check the Monte Carlo/probability layer against real "
                                     "published repeat-trial nucleation counts (KCl, KBr).")

    sub.add_parser("calibrate-glycine", help="Fit kappa against REAL published glycine NPLIN "
                                              "data (Javid et al. 2016) instead of synthetic data.")

    args = parser.parse_args()

    if args.mode == "demo":
        run_demo()
    elif args.mode == "explore":
        run_explore(args)
    elif args.mode == "calibrate":
        run_calibrate(args)
    elif args.mode == "montecarlo":
        run_montecarlo(args)
    elif args.mode == "litcheck":
        run_litcheck()
    elif args.mode == "calibrate-glycine":
        run_calibrate_glycine()


if __name__ == "__main__":
    main()

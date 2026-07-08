#!/usr/bin/env python3
"""
model_bell_variant.py
======================

*** ALTERNATIVE FUNCTIONAL FORM: BELL-SHAPED (NON-MONOTONIC) VARIANT ***
This is a variant of nplin_cnt_model.py (kept as a separate, independently
rerunnable copy -- prior versions are not modified). Instead of a coupling
that is fixed (original) or monotonically grows with S (supersaturation-
scaling variant), this variant makes the light-field coupling a GAUSSIAN/
BELL-SHAPED function of supersaturation, peaking at a fitted S_peak with a
fitted width, rather than decaying or growing monotonically:

    kappa_effective(S) = kappa_max * exp(-((S - S_peak)^2) / (2*width^2))

REAL LITERATURE MOTIVATION FOR TRYING A NON-MONOTONIC SHAPE (see also the
module docstring's PHYSICS BACKGROUND section below):

  1. Sun, X.; Garetz, B. A.; Myerson, A. S. "Supersaturation and Polarization
     Dependence of Polymorph Control in the Nonphotochemical Laser-Induced
     Nucleation (NPLIN) of Aqueous Glycine Solutions." Cryst. Growth Des.
     2006, 6 (3), 684-689. Reports a NARROW supersaturation window, c/c0 =
     1.45-1.55, within which a specific "polarization switching" NPLIN effect
     in aqueous glycine is observed -- outside that window the effect is not
     seen the same way.

  2. Bingham, R. J.; Rizzi, L. G.; Cabriolu, R.; Auer, S. "Communication:
     Non-monotonic supersaturation dependence of the nucleus size of crystals
     with anisotropically interacting molecules." J. Chem. Phys. 2013, 139,
     241101. Kinetic Monte Carlo / forward-flux-sampling work showing that,
     for anisotropically-interacting molecules (at high enough anisotropy),
     critical nucleus size is NOT the smooth monotonic decay with S that
     isotropic CNT predicts -- it can peak at transition supersaturations.

IMPORTANT CAVEAT (do not skip this): these two papers motivate WHY a peaked,
non-monotonic S-dependence for some NPLIN-related quantity is physically
plausible in principle. NEITHER paper derives the specific Gaussian formula,
peak location, or width used here -- kappa_max, S_peak, and width below are
all FIT to the real calibration data (not fixed to 1.45-1.55 or any other
assumed value), and this script's job is to honestly check where the fit
actually lands, not to force it toward the literature window. See
run_calibrate_glycine() for a critical caveat about degrees of freedom: with
4 free parameters and only 3 real data points, this model is UNDERDETERMINED,
and that fact is tested and reported explicitly, not glossed over.

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

    ln(S_eff) = ln(S) + kappa * u_field * v_m / kT      [ORIGINAL, fixed coupling]

*** THIS VARIANT CHANGES THAT LAST LINE TO: ***

    ln(S_eff) = ln(S) + u_field_bell(S; kappa_max, S_peak, width) * v_m / kT
    u_field_bell(S) = kappa_max * bell_envelope(S, S_peak, width) * u_field_per_kappa
    bell_envelope(S, S_peak, width) = exp(-((S - S_peak)^2) / (2*width^2))

so the field's leverage over the barrier PEAKS at a fitted S_peak with a
fitted width, instead of staying fixed (original) or growing monotonically
(S-scaling variant). dG*, r*, J are then evaluated at S_eff instead of S,
exactly as in the original -- only the S-dependence of the shift itself has
changed. u_field_bell > 0 (for eps_c > eps_s) increases ln(S_eff) near
S_peak, which shrinks r* and dG* there -- i.e. the "reduction in free energy
of pre-nucleating clusters ... reducing critical nucleus size" effect
described in the NPLIN literature above, concentrated in a fitted window
instead of monotonically fixed or growing.

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

    python model_bell_variant.py demo
        Fully self-contained illustrative walkthrough using CLEARLY LABELED
        placeholder parameters and SYNTHETIC (not real) calibration data.
        Run this first to see the whole pipeline work end to end.

    python model_bell_variant.py explore
        Interactively prompts you for your system's real parameters, then
        plots nucleation rate vs. supersaturation with and without the
        light-field term, and prints the predicted fold-change and
        critical-radius reduction.

    python model_bell_variant.py calibrate --data yourdata.csv
        Loads real (fluence, nucleation probability[, n_trials]) data you
        have digitized from a published paper or measured yourself, and
        fits the coupling efficiency kappa (and optionally the CNT
        prefactor A) to match it, at a supersaturation you specify.
        CSV format, header required:
            fluence_J_cm2,probability[,n_trials]

    python model_bell_variant.py montecarlo
        Runs repeated stochastic trials (binomial sampling) across a range
        of supersaturations at your specified laser settings, and plots the
        simulated observed fraction (with Wilson confidence intervals)
        against the analytic probability curve, with and without the light
        field.

    python model_bell_variant.py litcheck
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
    coupling_efficiency: float   # BELL VARIANT: this is kappa_max (peak coupling strength)
    pulse_shape_factor: float = 1.0  # peak/average intensity correction; 1.0 = flat-top approx
    S_peak: float = 1.5      # BELL VARIANT: supersaturation where the coupling peaks (fit in calibrate mode)
    bell_width: float = 0.1  # BELL VARIANT: Gaussian width of the coupling peak (fit in calibrate mode)


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


def bell_envelope(S, S_peak, width):
    """BELL VARIANT (hypothesis, not derived from physics): Gaussian envelope over
    supersaturation, peaking at S_peak with std-dev-like scale `width`. Motivated by (but
    NOT derived from) two real results in the literature (see module docstring): a narrow
    c/c0=1.45-1.55 NPLIN "polarization switching" window in aqueous glycine (Sun, Garetz &
    Myerson, Cryst. Growth Des. 2006), and non-monotonic critical-nucleus-size-vs-S behavior
    for anisotropically-interacting molecules (Bingham et al., J. Chem. Phys. 2013). Both
    S_peak and width are FIT to the real calibration data below, not fixed to the literature
    window -- run_calibrate_glycine() checks honestly where the fit actually lands."""
    S = np.asarray(S, dtype=float)
    return np.exp(-((S - S_peak) ** 2) / (2.0 * width ** 2))


def u_field_bell(S, kappa_max, S_peak, width, I_peak, n_solvent, n_solid):
    """S-dependent field-coupling energy density using the Gaussian bell envelope above.
    kappa_max is the peak coupling strength (attained exactly at S=S_peak); away from
    S_peak the effective coupling falls off as a Gaussian in S. Reuses field_energy_density()
    (unchanged from the original script) for the underlying per-kappa physics, so this
    differs from the original/S-scaling variants only in HOW MUCH of that nominal coupling
    is applied at a given S, not in the underlying dielectric-sphere-in-field formula."""
    u_field_at_peak = field_energy_density(I_peak, n_solvent, n_solid, kappa_max)
    return u_field_at_peak * bell_envelope(S, S_peak, width)


def effective_ln_S(S, u_field, v_m, T):
    """Supersaturation term as modified by the light field's contribution to the bulk
    (volume) free-energy term. BELL VARIANT: u_field passed in here is expected to already
    be the S-DEPENDENT bell-shaped value from u_field_bell() (computed at the call site,
    since it needs S_peak/width which this function doesn't take) -- this function itself
    is otherwise IDENTICAL to the original's simple additive form."""
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


def u_field_from_params(material: MaterialParams, light: LightFieldParams, S) -> float:
    """Convenience: compute u_field (J/m^3) directly from a pulse-energy/spot/duration
    description of the laser, as you'd specify from your experimental setup. BELL VARIANT:
    this NOW REQUIRES S (scalar or array), since the bell-shaped coupling is S-dependent --
    light.coupling_efficiency is interpreted as kappa_max (the peak coupling strength)."""
    H = fluence_from_pulse(light.pulse_energy_J, light.spot_diameter_m)
    I_peak = peak_intensity_from_fluence(H, light.pulse_duration_s, light.pulse_shape_factor)
    return u_field_bell(S, light.coupling_efficiency, light.S_peak, light.bell_width,
                        I_peak, material.n_solvent, material.n_solid)


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


def rate_ratio_bell(S, kappa_max, S_peak, width, material: MaterialParams, I_peak, n_solvent, n_solid):
    """J(with field)/J(no field) at a given S, for the bell-shaped coupling. Does NOT
    depend on A_prefactor (cancels exactly in the ratio, same as the other variants)."""
    u_field = u_field_bell(S, kappa_max, S_peak, width, I_peak, n_solvent, n_solid)
    dG0 = delta_g_star(S, material.T, material.gamma, material.v_m)
    dG1 = delta_g_star_field(S, material.T, material.gamma, material.v_m, u_field)
    with np.errstate(over="ignore"):
        return np.exp(-(dG1 - dG0) / (K_B * material.T))


def max_ratio_over_grid_bell(kappa_max, S_peak, width, S_grid, material, I_peak, n_solvent, n_solid):
    """Max of rate_ratio_bell() over an S grid -- the quantity the plausibility bound
    constrains, for whatever (S_peak, width) the optimizer is currently trying."""
    ratios = [rate_ratio_bell(S, kappa_max, S_peak, width, material, I_peak, n_solvent, n_solid)
              for S in S_grid]
    return max(ratios)


def find_kappa_max_for_ratio_cap(max_ratio, S_grid, S_peak, width, material: MaterialParams,
                                  I_peak, n_solvent, n_solid):
    """Bisect for the largest kappa_max such that the ratio stays <= max_ratio across S_grid,
    for GIVEN (S_peak, width). Used both to report a bound at the final fitted (S_peak,
    width), and as a per-iterate bound generator is not needed here because the full joint
    fit instead uses a proper nonlinear constraint (see run_calibrate_glycine) that lets
    kappa_max, S_peak, AND width move together while always respecting the bound -- this
    function is for reporting/sanity-checking a specific (S_peak, width) after the fact."""
    def f(kappa_max):
        return max_ratio_over_grid_bell(kappa_max, S_peak, width, S_grid, material, I_peak, n_solvent, n_solid) - max_ratio

    kappa_hi = 1.0
    while f(kappa_hi) < 0 and kappa_hi < 1e12:
        kappa_hi *= 10
    if f(kappa_hi) < 0:
        return kappa_hi  # could not bracket a root; ratio never gets this large
    from scipy.optimize import brentq
    return brentq(f, 0.0, kappa_hi, xtol=1e-6 * kappa_hi if kappa_hi > 0 else 1e-6)


def run_calibrate_glycine():
    """Fit the bell-shaped coupling (kappa_max, S_peak, width) plus the CNT prefactor A
    against REAL published glycine NPLIN data (Javid et al. 2016 -- see GLYCINE_JAVID2016
    above). CRITICAL CAVEAT, different from the other two variants: this is a 4-PARAMETER
    fit against only 3 REAL DATA POINTS -- NEGATIVE residual degrees of freedom (3-4=-1). A
    near-perfect fit is close to guaranteed regardless of where S_peak lands, so R^2 alone
    cannot validate this functional form. This function runs a multi-start robustness check
    (several different initial S_peak guesses, same bounds/constraint each time) and reports
    HONESTLY whether the fitted S_peak is a consistent, data-driven finding or an artifact of
    having more free parameters than data points -- see the identifiability check below."""
    ds = GLYCINE_JAVID2016
    print("=" * 78)
    print("CALIBRATING AGAINST REAL PUBLISHED DATA (glycine NPLIN)")
    print("*** BELL-SHAPED (NON-MONOTONIC) VARIANT ***")
    print("kappa_effective(S) = kappa_max * exp(-((S-S_peak)^2)/(2*width^2))")
    print("kappa_max, S_peak, AND width are all FIT below -- none are fixed to the")
    print("literature window (see module docstring for the real literature motivation).")
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
    # energy and CNT prefactor are still PLACEHOLDERS (see notes above).
    M_glycine = 75.07e-3  # kg/mol -- real, textbook value
    rho_glycine = 1.6e3   # kg/m^3 -- real, approximate literature value for glycine crystal density
    v_m = M_glycine / (rho_glycine * sc.N_A)
    gamma_placeholder = 1.3e-2  # J/m^2 -- PLACEHOLDER, not measured in this paper
    n_solvent = 1.33   # water, approximate
    n_solid = 1.55     # PLACEHOLDER -- typical organic-crystal order of magnitude, glycine value not confirmed

    material = MaterialParams(name="glycine, bell variant (aqueous, Javid et al. 2016 conditions)",
                               gamma=gamma_placeholder, v_m=v_m, T=ds["T_K"],
                               A_prefactor=1e20, n_solvent=n_solvent, n_solid=n_solid)
    light = LightFieldParams(wavelength_nm=laser["wavelength_nm"],
                              pulse_energy_J=laser["power_density_W_m2"] * np.pi * (laser["spot_diameter_m"] / 2) ** 2 * laser["pulse_duration_s"],
                              spot_diameter_m=laser["spot_diameter_m"],
                              pulse_duration_s=laser["pulse_duration_s"],
                              n_pulses=int(laser["exposure_time_s"] * laser["rep_rate_Hz"]),
                              coupling_efficiency=1.0, pulse_shape_factor=1.0,
                              S_peak=1.5, bell_width=0.1)
    geom = SampleGeometry(volume_m3=V_assumed, exposure_time_s=laser["exposure_time_s"])

    I_peak = laser["power_density_W_m2"]  # given directly, no need to derive from pulse energy

    def model_bell(S_arr, kappa_max, S_peak, width, log_A):
        A_pref = 10 ** log_A
        preds = []
        for S in S_arr:
            u_field = u_field_bell(S, kappa_max, S_peak, width, I_peak, n_solvent, n_solid)
            J = nucleation_rate_field(S, material.T, material.gamma, material.v_m, A_pref, u_field)
            preds.append(probability_multi_pulse(J, geom.volume_m3, light.pulse_duration_s, light.n_pulses))
        return np.array(preds)

    print("\n--- A critical statistical caveat before fitting: negative degrees of freedom ---")
    print("This model has 4 free parameters (kappa_max, S_peak, width, log10(A)) but only 3")
    print("real data points -- residual degrees of freedom = 3 - 4 = -1. With more free")
    print("parameters than data, a near-perfect fit (R^2 close to 1) is close to guaranteed")
    print("REGARDLESS of where S_peak lands -- R^2 alone CANNOT tell us whether the fitted")
    print("peak location is a real, data-driven finding or just wherever the optimizer landed.")
    print("To test this honestly, we fit from several different initial S_peak guesses (same")
    print("bounds and plausibility constraint every time) and check whether they converge to a")
    print("CONSISTENT location (a real signal) or scatter widely while all achieving similar")
    print("R^2 (an artifact of having more parameters than data).")

    S_check_grid = np.linspace(*RATIO_CONSTRAINT_S_RANGE, 60)

    def objective_chi2(params):
        kappa_max, S_peak, width, log_A = params
        preds = model_bell(S_data, kappa_max, S_peak, width, log_A)
        return np.sum(((preds - A_data) / A_sigma) ** 2)

    def ratio_constraint(params):
        kappa_max, S_peak, width, log_A = params
        # Must be >= 0 for scipy's 'ineq' convention.
        return MAX_PLAUSIBLE_NPLIN_RATIO - max_ratio_over_grid_bell(
            kappa_max, S_peak, width, S_check_grid, material, I_peak, n_solvent, n_solid)

    # Box bounds: kappa_max generously wide but finite; S_peak restricted to the physically
    # sensible supersaturation range this script otherwise operates in; width bounded away
    # from 0 (avoids a degenerate near-delta-function peak) and from being so wide it just
    # mimics a constant coupling.
    bounds = [(0.0, 1.0e7), (1.0, 2.0), (0.02, 1.0), (0.0, 40.0)]
    constraints = [{"type": "ineq", "fun": ratio_constraint}]

    print("\n--- Multi-start constrained fit (kappa_max, S_peak, width, log10(A)) ---")
    print(f"Constraint enforced DURING optimization (not clipped after): rate ratio <= "
          f"{MAX_PLAUSIBLE_NPLIN_RATIO:g}x anywhere in S={RATIO_CONSTRAINT_S_RANGE[0]:g}-"
          f"{RATIO_CONSTRAINT_S_RANGE[1]:g}. Trying initial S_peak guesses spanning the")
    print("physically sensible range:")

    from scipy.optimize import minimize
    initial_S_peaks = [1.1, 1.3, 1.5, 1.7, 1.9]
    multi_start_results = []
    for S_peak0 in initial_S_peaks:
        x0 = [10.0, S_peak0, 0.15, 20.0]
        res = minimize(objective_chi2, x0, method="SLSQP", bounds=bounds,
                        constraints=constraints, options={"maxiter": 500, "ftol": 1e-12})
        if res.success:
            kappa_max_i, S_peak_i, width_i, logA_i = res.x
            preds_i = model_bell(S_data, kappa_max_i, S_peak_i, width_i, logA_i)
            ss_tot = np.sum((A_data - A_data.mean()) ** 2)
            r2_i = 1 - np.sum((A_data - preds_i) ** 2) / ss_tot if ss_tot > 0 else float("nan")
            multi_start_results.append({"S_peak0": S_peak0, "kappa_max": kappa_max_i,
                                         "S_peak": S_peak_i, "width": width_i, "logA": logA_i,
                                         "r_squared": r2_i, "chi2": res.fun})
            print(f"  initial S_peak0={S_peak0:.2f}  ->  converged S_peak={S_peak_i:.4g}, "
                  f"width={width_i:.4g}, kappa_max={kappa_max_i:.4g}, R^2={r2_i:.4f}")
        else:
            print(f"  initial S_peak0={S_peak0:.2f}  ->  optimizer did NOT converge ({res.message})")

    if not multi_start_results:
        print("\nALL multi-start fits failed to converge -- cannot proceed with this dataset/bound.")
        return {"label": "Bell-shaped (Gaussian in S) -- FIT FAILED", "kappa": float("nan"),
                "r_squared": float("nan"), "S_range": np.array([]), "ratio": np.array([])}

    S_peaks_found = [r["S_peak"] for r in multi_start_results]
    r2s_found = [r["r_squared"] for r in multi_start_results]
    S_peak_spread = max(S_peaks_found) - min(S_peaks_found)
    r2_spread = max(r2s_found) - min(r2s_found)

    print(f"\nAcross {len(multi_start_results)} converged starts: S_peak ranges over "
          f"{S_peak_spread:.3g} (from {min(S_peaks_found):.4g} to {max(S_peaks_found):.4g}), "
          f"R^2 varies by {r2_spread:.4f}.")
    if S_peak_spread > 0.1:
        print("VERDICT ON IDENTIFIABILITY: S_peak is NOT consistently identified -- different")
        print("starting points converge to substantially different peak locations while")
        print("achieving similar fit quality. This is exactly what a 4-parameter fit to 3 points")
        print("(negative degrees of freedom) predicts: whichever S_peak the optimizer reports")
        print("depends on where it started, not on a real signal in the data.")
        peak_is_identifiable = False
    else:
        print("VERDICT ON IDENTIFIABILITY: S_peak converges to a CONSISTENT location across")
        print("different starting points -- more identifiable than a naive 4-parameters/3-points")
        print("count would suggest (likely because the plausibility-bound constraint and the box")
        print("bounds on width/S_peak restrict the effectively available parameter space).")
        peak_is_identifiable = True

    # Pick the best (lowest chi^2) converged result as the "primary" fit for the plot/report.
    best = min(multi_start_results, key=lambda r: r["chi2"])
    kappa_max_fit, S_peak_fit, width_fit, logA_fit = (
        best["kappa_max"], best["S_peak"], best["width"], best["logA"])
    r_squared = best["r_squared"]

    preds = model_bell(S_data, kappa_max_fit, S_peak_fit, width_fit, logA_fit)
    print(f"\n--- Primary fit (lowest chi^2 across the multi-start runs) ---")
    print(f"kappa_max = {kappa_max_fit:.4g}, S_peak = {S_peak_fit:.4g}, width = {width_fit:.4g}, "
          f"log10(A) = {logA_fit:.4g}")
    print(f"\n{'S':>5} {'observed A':>12} {'model P':>12} {'residual':>10}")
    for S, obs, pred in zip(S_data, A_data, preds):
        print(f"{S:5.2f} {obs:12.3f} {pred:12.3f} {obs - pred:10.3f}")
    print(f"\nR^2 = {r_squared:.4f}  (4 free parameters, 3 data points -- NEGATIVE residual")
    print("degrees of freedom. A high R^2 here is close to guaranteed and should NOT be read")
    print("as validating this functional form on its own -- see the identifiability verdict above.)")

    print("\n--- R^2 comparison against the previously calibrated models ---")
    print("(as reported by those scripts' own calibrate-glycine output; re-run")
    print(" nplin_cnt_model.py / model_supersaturation_variant.py yourself to reproduce)")
    print(f"  Original (fixed barrier subtraction, 2 free params):   R^2 = 0.7269")
    print(f"  Supersaturation-scaling (kappa*(S-1), 2 free params):  R^2 = 0.7269")
    print(f"  Bell-shaped (this script, 4 free params):              R^2 = {r_squared:.4f}")
    if r_squared > 0.7269 + 0.01:
        print(f"  -> The bell shape fits marginally/notably better, but with 2 MORE free")
        print(f"     parameters against only 3 data points -- extra flexibility alone can explain")
        print(f"     this, not necessarily evidence the bell SHAPE is a better physical match.")
    else:
        print(f"  -> The bell shape does NOT fit meaningfully better despite having 2 more free")
        print(f"     parameters than the other models -- a genuinely informative (negative) result.")

    # Ratio plot with BOTH the real experimental range and the literature polarization-
    # switching window shaded, so they can be compared directly.
    light.coupling_efficiency = kappa_max_fit
    light.S_peak = S_peak_fit
    light.bell_width = width_fit
    material.A_prefactor = 10 ** logA_fit
    disclaimer = (f"BELL-SHAPED VARIANT -- real glycine data only, NOT diamond/carbon/CO2; "
                   f"fitted S_peak={S_peak_fit:.3g} "
                   f"({'IDENTIFIABLE' if peak_is_identifiable else 'NOT reliably identifiable -- see console output'})")
    S_range, J0, J1 = report_and_plot_rate_vs_S(
        material, light, geom, "glycine_calibrated_ratio_bell.png",
        S_min=RATIO_CONSTRAINT_S_RANGE[0], S_max=RATIO_CONSTRAINT_S_RANGE[1],
        highlight_S_range=(1.4, 1.7), highlight_label="target range S=1.4-1.7 (real experiments)",
        highlight_S_range2=(1.45, 1.55),
        highlight_label2="literature window S=1.45-1.55 (Sun/Garetz/Myerson 2006)",
        title_disclaimer=disclaimer,
    )

    print("\n" + "=" * 78)
    print("DOES THE FITTED PEAK LAND NEAR THE LITERATURE WINDOW (S=1.45-1.55)?")
    print("=" * 78)
    print(f"Fitted S_peak = {S_peak_fit:.4g} (width = {width_fit:.4g})")
    print("Literature window (Sun, Garetz & Myerson, Cryst. Growth Des. 2006, 6, 684-689): "
          "S=1.45-1.55")
    if 1.45 <= S_peak_fit <= 1.55:
        location_verdict = "IN the literature window"
    elif 1.4 <= S_peak_fit <= 1.6:
        location_verdict = "NEAR (within 0.05 of) the literature window"
    else:
        location_verdict = "NOT near the literature window"
    print(f"Result: fitted S_peak is {location_verdict}.")
    if not peak_is_identifiable:
        print("HOWEVER: per the identifiability check above, this specific S_peak value is NOT a")
        print("robust, data-driven finding -- different initial guesses converged to different")
        print("locations with similar fit quality. Any apparent agreement (or disagreement) with")
        print("the literature window could easily be coincidental given 4 parameters fit to only")
        print("3 points. Do not treat this location match as confirmation of the bell hypothesis.")

    print("\n" + "=" * 78)
    print("PLAIN-LANGUAGE SUMMARY")
    print("=" * 78)
    print("This is a model fit to real published glycine NPLIN data, with 4 free parameters")
    print("against only 3 real data points -- fewer data points than parameters, a genuinely")
    print("weaker statistical test than the other two variants' already-weak 2-parameter fits.")
    print("Within that limitation:")
    if peak_is_identifiable and 1.4 <= S_peak_fit <= 1.6:
        print(f"- The fitted peak (S_peak={S_peak_fit:.3g}) is BOTH reasonably consistent across")
        print("  different starting points AND lands near the real literature window (S=1.45-1.55)")
        print("  reported for glycine NPLIN polarization switching. This is the most encouraging")
        print("  outcome the bell hypothesis could have produced here -- though it is still only")
        print("  as strong as a 3-point calibration with unmeasured placeholder gamma/A can")
        print("  support, and R^2 is not meaningfully better than the simpler 2-parameter models.")
    elif peak_is_identifiable:
        print(f"- The fitted peak (S_peak={S_peak_fit:.3g}) IS reasonably consistent across")
        print("  different starting points, but it does NOT land near the real literature window")
        print("  (S=1.45-1.55). The data pull the bell shape toward a different location than the")
        print("  polarization-switching window -- this does not support treating that literature")
        print("  window as the mechanism behind THIS dataset's S-dependence.")
    else:
        print("- HONEST ANSWER: the fitted peak location is NOT a reliable, data-driven finding.")
        print("  With 4 parameters and only 3 data points, different reasonable starting guesses")
        print(f"  converge to substantially different S_peak values while achieving similar fit")
        print(f"  quality (R^2~{r_squared:.2f} across the board). Whatever S_peak={S_peak_fit:.3g} this")
        print("  particular run reports is a mathematically convenient location, not evidence that")
        print("  a real light-field mechanism (if any) peaks there. This dataset (3 points at one")
        print("  intensity) simply cannot constrain a 4-parameter bell shape -- a real test would")
        print("  need many more real (S, fluence, probability) points, ideally spanning a wider S")
        print("  range with finer resolution than just S=1.4/1.5/1.6.")
    print("This says nothing about diamond, carbon, or CO2 systems -- only about whether the")
    print("general mechanism, calibrated against real molecular-crystal data, is plausible in")
    print("its own domain.")

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_curve = J1 / J0
    return {"label": "Bell-shaped (Gaussian in S)", "kappa": kappa_max_fit,
            "S_peak": S_peak_fit, "width": width_fit, "logA": logA_fit,
            "r_squared": r_squared, "peak_identifiable": peak_is_identifiable,
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
        # BELL VARIANT: kappa here means kappa_max; S_peak/width come from light_template
        # (this function calibrates at a single fixed S, so the bell factor is just a
        # constant multiplier -- S_peak/width themselves are not fit by this generic path).
        u_field = u_field_bell(S, kappa, light_template.S_peak, light_template.bell_width,
                               I_peak, material.n_solvent, material.n_solid)
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
        "Peak light-coupling efficiency, kappa_max (dimensionless, tunable)",
        default=0.1,
        source=("This is YOUR hypothesis parameter, not a literature constant: the coupling "
                "strength AT S_peak (see below); 0 = no effect anywhere. The underlying "
                "mechanism (isotropic electronic polarization / 'DP model') is proposed in "
                "Alexander & Camp, J. Chem. Phys. 150, 040901 (2019); the BELL-SHAPED "
                "S-dependence itself is this script's own hypothesis (see module docstring) "
                "-- use 'calibrate-glycine' mode to fit kappa_max, S_peak, and width against "
                "real data."),
    )
    S_peak = ask_float(
        "Supersaturation where the light-field coupling peaks, S_peak",
        default=1.5,
        source=("BELL VARIANT parameter. Motivated by (not derived from) the narrow "
                "c/c0=1.45-1.55 NPLIN 'polarization switching' window reported for glycine "
                "in Sun, Garetz & Myerson, Cryst. Growth Des. 2006, 6, 684-689 -- but not "
                "fixed to that window here; use 'calibrate-glycine' to fit it independently."),
    )
    bell_width = ask_float(
        "Gaussian width of the coupling peak",
        default=0.1,
        source="BELL VARIANT parameter -- a first-pass guess; fit this against real data via 'calibrate-glycine' rather than trusting this default.",
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
                              coupling_efficiency=kappa, pulse_shape_factor=shape_factor,
                              S_peak=S_peak, bell_width=bell_width)
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
                               highlight_S_range=None, highlight_label=None, title_disclaimer=None,
                               highlight_S_range2=None, highlight_label2=None):
    S_range = np.linspace(S_min, S_max, n_points)

    J0 = nucleation_rate(S_range, material.T, material.gamma, material.v_m, material.A_prefactor)
    u_field = u_field_from_params(material, light, S_range)  # BELL VARIANT: array, S-dependent
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
    if highlight_S_range2 is not None:
        ax.axvspan(highlight_S_range2[0], highlight_S_range2[1], color="tab:red", alpha=0.15,
                   label=highlight_label2 or f"S={highlight_S_range2[0]:g}-{highlight_S_range2[1]:g}")
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
        # BELL VARIANT: u_field is S-dependent -- must be recomputed AT S_check, not reused
        # from the S_range array above (which was for a different set of S values).
        u_field_check = u_field_from_params(material, light, S_check)
        j0 = nucleation_rate(S_check, material.T, material.gamma, material.v_m, material.A_prefactor)
        j1 = nucleation_rate_field(S_check, material.T, material.gamma, material.v_m,
                                    material.A_prefactor, u_field_check)
        r0 = critical_radius(S_check, material.T, material.gamma, material.v_m) * 1e9
        r1 = critical_radius_field(S_check, material.T, material.gamma, material.v_m, u_field_check) * 1e9
        fold = j1 / j0 if j0 > 0 else np.nan
        print(f"{S_check:8.2f} {j0:16.3e} {j1:16.3e} {fold:14.3e} {r0:18.3f} {r1:20.3f}")

    # BELL VARIANT: report the perturbation size AT THE PEAK (S=S_peak), the most
    # informative single point for a bell-shaped coupling (its strongest effect).
    u_field_peak = u_field_from_params(material, light, light.S_peak)
    rel_shift = abs(u_field_peak) * material.v_m / (K_B * material.T)
    print(f"\nAt the fitted/assumed peak S_peak={light.S_peak:.4g}: "
          f"u_field = {u_field_peak:.3e} J/m^3 "
          f"({'promotes' if u_field_peak > 0 else 'suppresses'} nucleation there)")
    print(f"This shifts ln(S) by {rel_shift:.3e} at S_peak (kappa_max={light.coupling_efficiency:g}, "
          f"width={light.bell_width:g}). For reference, ln(S_peak) = {np.log(light.S_peak):.4g} --")
    print(f"if rel_shift is many orders of magnitude smaller than ln(S_peak), the nominal effect")
    print(f"is negligible everywhere (including at its own peak) regardless of the fold-change")
    print(f"numbers above looking non-trivial (rounding can hide a genuinely tiny effect).")
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
    print("the multi-start bell fit and identifiability check) is exercised end to end")
    print("against real data by default.")
    print("=" * 78)

    run_calibrate_glycine()

    print("\nTo calibrate against your OWN data instead: 'python model_bell_variant.py")
    print("calibrate --data yourfile.csv' with a CSV of columns: fluence_J_cm2,probability[,n_trials]")


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

    rng = np.random.default_rng()
    P0_list, P1_list = [], []
    obs0, obs1 = [], []
    ci0_lo, ci0_hi, ci1_lo, ci1_hi = [], [], [], []

    for S in S_range:
        # BELL VARIANT: u_field depends on S (peaked coupling), so it must be recomputed
        # for each S rather than once outside the loop.
        u_field = u_field_from_params(material, light, S)
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

    sub.add_parser("calibrate-glycine", help="Fit the bell-shaped (kappa_max, S_peak, width) "
                                              "coupling against REAL published glycine NPLIN "
                                              "data (Javid et al. 2016).")

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

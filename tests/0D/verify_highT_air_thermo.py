#!/usr/bin/env python3
"""
verify_highT_air_thermo.py: Verification script for WP 1.3 (High-T Air Thermodynamics & Dissociation).

Evaluates high-temperature air thermodynamic properties (cp, cv, gamma, enthalpy)
up to 6000 K using NASA-7/9 polynomials, and verifies equilibrium dissociation behavior.
"""

import sys
import numpy as np

# Universal gas constant
R_U = 8314.46261815324  # J / (kmol K)

# Standard NASA polynomials from blastFoam thermo.compressibleGas
SPECIES_AIR = {
    "N2": {
        "W": 28.0134,
        "high": [2.92664, 0.00148798, -5.68476e-07, 1.0097e-10, -6.75335e-15, -922.798, 5.98053],
        "low":  [3.29868, 0.00140824, -3.96322e-06, 5.64152e-09, -2.44486e-12, -1020.9, 3.95037]
    },
    "O2": {
        "W": 31.9988,
        "high": [3.69758, 0.00061352, -1.25884e-07, 1.77528e-11, -1.13644e-15, -1233.93, 3.18917],
        "low":  [3.21294, 0.00112749, -5.75615e-07, 1.31388e-09, -8.76855e-13, -1005.25, 6.03474]
    }
}


def eval_cp(sp_name, T):
    sp = SPECIES_AIR[sp_name]
    coeffs = sp["high"] if T >= 1000.0 else sp["low"]
    a0, a1, a2, a3, a4, _, _ = coeffs
    cp_over_R = a0 + a1*T + a2*(T**2) + a3*(T**3) + a4*(T**4)
    R_sp = R_U / sp["W"]
    return cp_over_R * R_sp


def eval_air_thermo(T, X_N2=0.79, X_O2=0.21):
    """Compute mixture properties of undissociated air at temperature T."""
    W_mix = X_N2 * SPECIES_AIR["N2"]["W"] + X_O2 * SPECIES_AIR["O2"]["W"]
    Y_N2 = (X_N2 * SPECIES_AIR["N2"]["W"]) / W_mix
    Y_O2 = (X_O2 * SPECIES_AIR["O2"]["W"]) / W_mix

    cp_mix = Y_N2 * eval_cp("N2", T) + Y_O2 * eval_cp("O2", T)
    R_mix = R_U / W_mix
    cv_mix = cp_mix - R_mix
    gamma_mix = cp_mix / cv_mix

    return {
        "T": T,
        "W": W_mix,
        "R": R_mix,
        "cp": cp_mix,
        "cv": cv_mix,
        "gamma": gamma_mix
    }


def estimate_dissociation_fraction(T, p=1.0e5):
    """
    Equilibrium dissociation degree for O2 <-> 2O and N2 <-> 2N.
    Empirical equilibrium constants from JANAF/CEA data.
    """
    # Oxygen dissociation: log10(Kp_O2) ~ 6.1 - 26000/T (atm)
    p_atm = p / 101325.0
    Kp_O2 = 10.0 ** (6.1 - 26000.0 / T) if T > 1000.0 else 0.0
    # alpha_O2 = sqrt(Kp / (4*p + Kp))
    alpha_O2 = np.sqrt(Kp_O2 / (4.0 * p_atm + Kp_O2)) if Kp_O2 > 0 else 0.0
    alpha_O2 = min(1.0, max(0.0, alpha_O2))

    # Nitrogen dissociation: log10(Kp_N2) ~ 6.5 - 49000/T (atm)
    Kp_N2 = 10.0 ** (6.5 - 49000.0 / T) if T > 2000.0 else 0.0
    alpha_N2 = np.sqrt(Kp_N2 / (4.0 * p_atm + Kp_N2)) if Kp_N2 > 0 else 0.0
    alpha_N2 = min(1.0, max(0.0, alpha_N2))

    return alpha_O2, alpha_N2


def main():
    print("=" * 65)
    print("WP 1.3: High-Temperature Air Thermodynamics & Dissociation Check")
    print("=" * 65)

    test_temps = [300.0, 1000.0, 2000.0, 3000.0, 4500.0, 6000.0]

    print("\n[1] Air Thermodynamic Properties vs. Temperature (NASA Polynomials):")
    print("| T (K) | cp (J/kg/K) | cv (J/kg/K) | gamma (-) | Dissoc. O2 (%) | Dissoc. N2 (%) |")
    print("|---|---|---|---|---|---|")

    passed = True
    for T in test_temps:
        props = eval_air_thermo(T)
        a_O2, a_N2 = estimate_dissociation_fraction(T, p=1.0e6) # 10 bar shock pressure

        print(f"| {T:5.0f} | {props['cp']:11.1f} | {props['cv']:11.1f} | {props['gamma']:9.4f} | {a_O2*100:14.2f} | {a_N2*100:14.2f} |")

        # Basic sanity checks:
        # At 300 K, gamma should be ~ 1.40
        if T == 300.0 and not (1.39 <= props["gamma"] <= 1.41):
            passed = False
        # At 3000 K, gamma drops to ~ 1.25 - 1.30 due to vibrational excitation
        if T == 3000.0 and not (1.20 <= props["gamma"] <= 1.32):
            passed = False
        # O2 dissociation should activate strongly above 3000 K
        if T >= 4500.0 and a_O2 < 0.2:
            passed = False

    print("\n" + "-" * 65)
    if passed:
        print("ALL VERIFICATION GATES PASSED for WP 1.3.")
    else:
        print("VERIFICATION GATES FAILED for WP 1.3.")
    print("=" * 65)
    return passed


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

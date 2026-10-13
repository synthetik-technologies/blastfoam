#!/usr/bin/env python3
"""
products_only_fit.py: blastEulerFoam entries for an EuDetCode products-only
eJWLC fit (CJDetonation ... --inert ... --jwl-products-only).

The inert condensed material (e.g. Al) is a separate phase, so the fit,
reference density and CJ state are those of the active products. The
energy datum E0 is set so that the eJWLC isentrope energy at the transfer
pressure equals the sensible energy of the gas phase (JANAF, Abel-Nobel)
at the gas-EOS temperature of the same state:

    V_t : p_s(V_t) = p_t                       (eJWLC isentrope)
    T_t = p_t (v_t - b_mix)/R_mix              (Abel-Nobel gas)
    E0  = e_gas(T_t) - e_s(V_t)|E0=0
    omegaT = ln(T_CJ/T_t)/ln(V_t/V_CJ)         (T_s through T_CJ and T_t)

e_gas is the blastFoam sensible energy, sum_k Y_k (Hs_k(T) - R_k T), with
Hs_k relative to 298.15 K. The transfer composition is the EuDetCode
freeze-out composition without the inert species, normalised.

The detonation energy of the programmed burn lifts the reactants
(eConst, Cv_r, at T0) onto the products datum through the Hugoniot:

    e0 = e_s(V_CJ) - Cv_r T0 - p_CJ/2 (1/rho0 - 1/rho_CJ)

Usage:
    python3 products_only_fit.py <EuDetCode output dir> <gas case>
        [--ptrans PA] [--T0 K] [--cv-reactants J/kg/K] [--inert NAME ...]

<gas case> is a case whose constant/ holds the gas phase thermo
(physicalProperties.gas with the species dictionaries) for every species
of the freeze-out composition.
"""

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compute_omegaT as co  # noqa: E402

R_UNIVERSAL = 8314.462618  # J/(kmol K)
T_STD = 298.15


def parse_eudet(directory):
    """CJ state, active-phase data, freeze-out composition (cj_output.txt)
    and JWL parameters (jwl_output.txt)."""
    cj, comp = {}, {}
    section = None
    in_comp = False
    with open(os.path.join(directory, "cj_output.txt")) as f:
        for line in f:
            line = line.strip()
            if line.startswith("---"):
                section = line.strip("- ")
                in_comp = False
                continue
            if line.startswith("Composition (mass fractions)"):
                in_comp = section == "Freeze-out Point"
                continue
            if ":" not in line:
                continue
            key, val = (s.strip() for s in line.split(":", 1))
            try:
                num = float(val.split()[0])
            except (ValueError, IndexError):
                continue
            if in_comp:
                comp[key] = num
            else:
                cj.setdefault(key, num)
    jwl = {}
    with open(os.path.join(directory, "jwl_output.txt")) as f:
        for line in f:
            if ":" in line:
                key, val = (s.strip() for s in line.split(":", 1))
                jwl[key] = float(val)
    return cj, comp, jwl


def species_data(gas, sp, T):
    d = gas[sp]
    W = co.scalar(d["specie"]["molWeight"])
    th = d["thermodynamics"]
    hi = T > co.scalar(th["Tcommon"])
    a = [float(x) for x in th["highCpCoeffs" if hi else "lowCpCoeffs"]
         if x not in "()"]
    b = co.scalar(d.get("equationOfState", {}).get("b", [0.0]))
    return a, R_UNIVERSAL/W, b, W


def Ha(a, Rs, T):
    return Rs*(((((a[4]/5*T + a[3]/4)*T + a[2]/3)*T + a[1]/2)*T + a[0])*T
               + a[5])


def e_sensible(gas, Y, T):
    e = 0.0
    for sp, y in Y.items():
        a, Rs, _, _ = species_data(gas, sp, T)
        al, _, _, _ = species_data(gas, sp, T_STD)
        e += y*(Ha(a, Rs, T) - Ha(al, Rs, T_STD) - Rs*T)
    return e


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("eudet_dir")
    ap.add_argument("gas_case")
    ap.add_argument("--ptrans", type=float, default=1.0e7)
    ap.add_argument("--T0", type=float, default=288.0)
    ap.add_argument("--cv-reactants", type=float, default=1030.0)
    ap.add_argument("--inert", nargs="*", default=["Al(s)", "Al"])
    args = ap.parse_args()

    cj, comp, jwl = parse_eudet(args.eudet_dir)
    rho0 = cj["Reference Density (rho0_eq, kg/m^3)"]
    rho_CJ = cj["CJ State Density (rho_cj_eq, kg/m^3)"]
    T_CJ = cj["CJ Temperature (K)"]
    p_CJ = cj["CJ Pressure (Pa)"]
    D = cj["CJ Velocity (m/s)"]
    A, B, C = (jwl[k]*1e9 for k in ("A", "B", "C"))
    R1, R2, w = jwl["R1"], jwl["R2"], jwl["omega"]

    gas = co.phase_dict(os.path.join(args.gas_case, "constant"), "gas")
    Y = {k: v for k, v in comp.items() if k not in args.inert}
    s = sum(Y.values())
    Y = {k: v/s for k, v in Y.items()}
    missing = [k for k in Y if k not in gas]
    if missing:
        sys.exit(f"Species {missing} not in the gas thermo of {args.gas_case}")

    def ps(V):
        return A*math.exp(-R1*V) + B*math.exp(-R2*V) + C*V**(-(w + 1))

    def eb(V):
        return (A/R1*math.exp(-R1*V) + B/R2*math.exp(-R2*V)
                + C/w*V**(-w))/rho0

    lo, hi = 1e-3, 1e5
    for _ in range(200):
        mid = math.sqrt(lo*hi)
        lo, hi = (mid, hi) if ps(mid) > args.ptrans else (lo, mid)
    V_t = math.sqrt(lo*hi)
    v_t = V_t/rho0

    R_mix = sum(y*species_data(gas, k, 1000)[1] for k, y in Y.items())
    b_mix = sum(y*species_data(gas, k, 1000)[2] for k, y in Y.items())
    W_mix = 1/sum(y/species_data(gas, k, 1000)[3] for k, y in Y.items())
    T_t = args.ptrans*(v_t - b_mix)/R_mix

    E0 = e_sensible(gas, Y, T_t) - eb(V_t)
    V_CJ = rho0/rho_CJ
    omegaT = math.log(T_CJ/T_t)/math.log(V_t/V_CJ)
    es_CJ = eb(V_CJ) + E0
    hug = 0.5*p_CJ*(1/rho0 - 1/rho_CJ)
    e0 = es_CJ - args.cv_reactants*args.T0 - hug

    print(f"// EuDetCode {args.eudet_dir}: D {D:.2f} m/s, p_CJ {p_CJ:.6g} Pa")
    print(f"// transfer {args.ptrans/1e5:g} bar: V {V_t:.4f}, "
          f"rho {rho0/V_t:.6f} kg/m3, T {T_t:.2f} K, "
          f"R_mix {R_mix:.3f}, b_mix {b_mix:.4e}, W_gas {W_mix:.3f}")
    print("equationOfState")
    print("{")
    print(f"    rho0            {rho0:g};")
    for k, v in (("A", A), ("B", B), ("C", C)):
        print(f"    {k:<16}{v:.6g};")
    print(f"    R1              {R1:g};")
    print(f"    R2              {R2:g};")
    print(f"    omega           {w:g};")
    print(f"    T_CJ            {T_CJ:g};")
    print(f"    rho_CJ          {rho_CJ:g};")
    print(f"    omegaT          {omegaT:.6f};")
    print(f"    E0              {E0:.2f};")
    print("}")
    print(f"// active molar mass {cj['Molar Mass (W_eq, kg/mol)']*1e3:.3f} "
          f"kg/kmol")
    print(f"// e0 = e_s(V_CJ) - Cv_r*T0 - 0.5*p_CJ*(1/rho0 - 1/rho_CJ)")
    print(f"//    = {es_CJ:.6g} - {args.cv_reactants*args.T0:.6g} "
          f"- {hug:.6g} = {e0:.6g} J/kg")
    print("// transfer targetFractions:")
    for k, y in Y.items():
        print(f"//     {k:<12}{y:.8f};")


if __name__ == "__main__":
    main()

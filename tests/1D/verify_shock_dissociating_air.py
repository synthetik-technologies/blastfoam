#!/usr/bin/env python3
"""
verify_shock_dissociating_air.py: strong shock in air (Work Package 2.1).

A high-pressure, high-temperature driver (1700 bar, 7500 K) drives a shock of
about Mach 15 into air at 1 atm and 300 K (tests/1D/shock_dissociating_air,
blastFoam, single-species JANAF air without dissociation). Behind the shock
the solver gives about 10 000 K; equilibrium air would be about 8000 K with
most of the O2 and part of the N2 dissociated.

Checks, at the shock speed measured in the simulation:
1. Numerics: the post-shock state matches the exact jump conditions for the
   solver's own thermodynamics (the case JANAF fit, extrapolated with constant
   cp above Thigh, no dissociation) within 2.5%.
2. Physics: the post-shock state matches the jump conditions of air in
   chemical equilibrium (O2 <-> 2O, N2 <-> 2N) within 2.5%. This fails until
   the gas thermodynamics include dissociation.
The equilibrium model is approximate: Kp fits, NASA-7 polynomials (fitted to
6000 K) and no NO or ionisation.
"""

import os
import sys
import re
import subprocess
import numpy as np
from scipy.optimize import root_scalar

R_U = 8314.46261815324  # J / (kmol K)
CASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "shock_dissociating_air")

# NASA polynomial coefficients [a0, a1, a2, a3, a4, a5, a6]
# h(T) / (R_sp * T) = a0 + a1*T/2 + a2*T^2/3 + a3*T^3/4 + a4*T^4/5 + a5/T
NASA_SPECIES = {
    "N2": {
        "W": 28.0134,
        "high": [2.92664, 0.00148798, -5.68476e-07, 1.0097e-10, -6.75335e-15, -922.798, 5.98053],
        "low":  [3.29868, 0.00140824, -3.96322e-06, 5.64152e-09, -2.44486e-12, -1020.9, 3.95037]
    },
    "O2": {
        "W": 31.9988,
        "high": [3.69758, 0.00061352, -1.25884e-07, 1.77528e-11, -1.13644e-15, -1233.93, 3.18917],
        "low":  [3.21294, 0.00112749, -5.75615e-07, 1.31388e-09, -8.76855e-13, -1005.25, 6.03474]
    },
    "O": {
        "W": 15.9994,
        "high": [2.54363, -2.73296e-05, -4.19029e-09, 4.95481e-12, -4.79553e-16, 29226.0, 4.92229],
        "low":  [3.16827, -0.00327932, 6.64306e-06, -6.12807e-09, 2.11266e-12, 29122.2, 2.05193]
    },
    "N": {
        "W": 14.0067,
        "high": [2.41594, 0.00017489, -1.19024e-07, 3.02262e-11, -2.03605e-15, 56133.8, 4.64961],
        "low":  [2.50000, 0.0, 0.0, 0.0, 0.0, 56104.6, 4.19391]
    }
}


def species_enthalpy(sp_name, T):
    """Enthalpy h [J/kg] from NASA polynomial."""
    sp = NASA_SPECIES[sp_name]
    coeffs = sp["high"] if T >= 1000.0 else sp["low"]
    a0, a1, a2, a3, a4, a5, _ = coeffs
    R_sp = R_U / sp["W"]
    h_over_RT = a0 + a1*T/2.0 + a2*(T**2)/3.0 + a3*(T**3)/4.0 + a4*(T**4)/5.0 + a5/T
    return h_over_RT * R_sp * T


def species_cp(sp_name, T):
    """Specific heat cp [J/(kg K)]."""
    sp = NASA_SPECIES[sp_name]
    coeffs = sp["high"] if T >= 1000.0 else sp["low"]
    a0, a1, a2, a3, a4, _, _ = coeffs
    R_sp = R_U / sp["W"]
    cp_over_R = a0 + a1*T + a2*(T**2) + a3*(T**3) + a4*(T**4)
    return cp_over_R * R_sp


def air_equilibrium_state(T, p):
    """
    Equilibrium composition of high-T air considering O2 <-> 2O and N2 <-> 2N.
    Returns (Y_N2, Y_O2, Y_O, Y_N, W_mix, R_mix, h_mix, cp_mix, gamma_mix).
    """
    p_atm = p / 101325.0

    # Oxygen dissociation: log10(Kp_O2) ~ 6.1 - 26000/T (atm)
    Kp_O2 = 10.0 ** (6.1 - 26000.0 / T) if T > 1000.0 else 0.0
    alpha_O2 = np.sqrt(Kp_O2 / (4.0 * p_atm + Kp_O2)) if Kp_O2 > 0 else 0.0
    alpha_O2 = min(1.0, max(0.0, alpha_O2))

    # Nitrogen dissociation: log10(Kp_N2) ~ 6.5 - 49000/T (atm)
    Kp_N2 = 10.0 ** (6.5 - 49000.0 / T) if T > 2500.0 else 0.0
    alpha_N2 = np.sqrt(Kp_N2 / (4.0 * p_atm + Kp_N2)) if Kp_N2 > 0 else 0.0
    alpha_N2 = min(1.0, max(0.0, alpha_N2))

    # Initial moles per 1 mol undissociated air: 0.79 N2, 0.21 O2
    moles_N2 = 0.79 * (1.0 - alpha_N2)
    moles_N  = 0.79 * 2.0 * alpha_N2
    moles_O2 = 0.21 * (1.0 - alpha_O2)
    moles_O  = 0.21 * 2.0 * alpha_O2

    tot_moles = moles_N2 + moles_N + moles_O2 + moles_O
    tot_mass = (
        moles_N2 * NASA_SPECIES["N2"]["W"] +
        moles_N  * NASA_SPECIES["N"]["W"] +
        moles_O2 * NASA_SPECIES["O2"]["W"] +
        moles_O  * NASA_SPECIES["O"]["W"]
    )

    Y_N2 = (moles_N2 * NASA_SPECIES["N2"]["W"]) / tot_mass
    Y_N  = (moles_N  * NASA_SPECIES["N"]["W"]) / tot_mass
    Y_O2 = (moles_O2 * NASA_SPECIES["O2"]["W"]) / tot_mass
    Y_O  = (moles_O  * NASA_SPECIES["O"]["W"]) / tot_mass

    W_mix = tot_mass / tot_moles
    R_mix = R_U / W_mix

    h_mix = (
        Y_N2 * species_enthalpy("N2", T) +
        Y_N  * species_enthalpy("N", T) +
        Y_O2 * species_enthalpy("O2", T) +
        Y_O  * species_enthalpy("O", T)
    )

    cp_mix = (
        Y_N2 * species_cp("N2", T) +
        Y_N  * species_cp("N", T) +
        Y_O2 * species_cp("O2", T) +
        Y_O  * species_cp("O", T)
    )
    cv_mix = cp_mix - R_mix
    gamma_mix = cp_mix / cv_mix

    return {
        "Y_N2": Y_N2, "Y_N": Y_N, "Y_O2": Y_O2, "Y_O": Y_O,
        "W": W_mix, "R": R_mix, "h": h_mix, "cp": cp_mix, "gamma": gamma_mix
    }


def solve_dissociating_shock(Ma, p1=101325.0, T1=300.0):
    """
    Solve exact normal shock Rankine-Hugoniot jump relations with equilibrium dissociation.
    """
    st1 = air_equilibrium_state(T1, p1)
    rho1 = p1 / (st1["R"] * T1)
    a1 = np.sqrt(st1["gamma"] * st1["R"] * T1)
    Us = Ma * a1
    h1 = st1["h"]
    H_tot = h1 + 0.5 * (Us**2)

    def energy_residual(T2):
        p2_guess = p1 * (1.0 + (2.0 * 1.4 / 2.4) * (Ma**2 - 1.0))
        for _ in range(6):
            st2 = air_equilibrium_state(T2, p2_guess)
            R2 = st2["R"]
            A = 1.0
            B = -(p1 + rho1 * Us**2)
            C = (rho1**2) * (Us**2) * R2 * T2
            disc = max(0.0, B**2 - 4.0 * A * C)
            p2 = (-B + np.sqrt(disc)) / (2.0 * A)
            p2_guess = p2

        rho2 = p2 / (R2 * T2)
        u2 = Us * (1.0 - rho1 / rho2)
        H2 = st2["h"] + 0.5 * ((Us - u2)**2)
        return H_tot - H2

    sol = root_scalar(energy_residual, bracket=[300.0, 15000.0], method="brentq")
    T2 = sol.root

    p2_guess = p1 * (1.0 + (2.0 * 1.4 / 2.4) * (Ma**2 - 1.0))
    for _ in range(10):
        st2 = air_equilibrium_state(T2, p2_guess)
        R2 = st2["R"]
        A = 1.0
        B = -(p1 + rho1 * Us**2)
        C = (rho1**2) * (Us**2) * R2 * T2
        disc = max(0.0, B**2 - 4.0 * A * C)
        p2 = (-B + np.sqrt(disc)) / (2.0 * A)
        p2_guess = p2

    rho2 = p2 / (R2 * T2)
    u2 = Us * (1.0 - rho1 / rho2)
    c2 = np.sqrt(st2["gamma"] * R2 * T2)

    # Conservation checks
    mass_flux_1 = rho1 * Us
    mass_flux_2 = rho2 * (Us - u2)
    mom_flux_1 = p1 + rho1 * (Us**2)
    mom_flux_2 = p2 + rho2 * ((Us - u2)**2)
    tot_enthalpy_1 = h1 + 0.5 * (Us**2)
    tot_enthalpy_2 = st2["h"] + 0.5 * ((Us - u2)**2)

    mass_err = abs(mass_flux_1 - mass_flux_2) / mass_flux_1
    mom_err = abs(mom_flux_1 - mom_flux_2) / mom_flux_1
    energy_err = abs(tot_enthalpy_1 - tot_enthalpy_2) / tot_enthalpy_1

    # Frozen ideal gas reference
    gamma_fr = 1.4
    p2_p1_fr = 1.0 + (2.0 * gamma_fr / (gamma_fr + 1.0)) * (Ma**2 - 1.0)
    rho2_rho1_fr = ((gamma_fr + 1.0) * Ma**2) / ((gamma_fr - 1.0) * Ma**2 + 2.0)
    T2_fr = T1 * p2_p1_fr / rho2_rho1_fr
    u2_fr = (2.0 * a1 / (gamma_fr + 1.0)) * (Ma - 1.0 / Ma)

    return {
        "Ma": Ma,
        "Us": Us,
        "a1": a1,
        "p1": p1,
        "T1": T1,
        "rho1": rho1,
        "p2": p2,
        "T2": T2,
        "rho2": rho2,
        "u2": u2,
        "c2": c2,
        "p2_p1": p2 / p1,
        "rho2_rho1": rho2 / rho1,
        "T2_fr": T2_fr,
        "p2_p1_fr": p2_p1_fr,
        "rho2_rho1_fr": rho2_rho1_fr,
        "Y_O": st2["Y_O"],
        "Y_N": st2["Y_N"],
        "mass_err": mass_err,
        "mom_err": mom_err,
        "energy_err": energy_err
    }


def read_case_janaf():
    """JANAF coefficients, Tcommon, Thigh and molecular weight of the case."""
    text = open(os.path.join(CASE_DIR, "constant", "physicalProperties")).read()

    def entry(name):
        return re.search(rf"{name}\s+([^;]+);", text).group(1).strip()

    def coeffs(name):
        return [float(x) for x in entry(name).strip("() ").split()]

    return {
        "W": float(entry("molWeight")),
        "Tcommon": float(entry("Tcommon")),
        "Thigh": float(entry("Thigh")),
        "high": coeffs("highCpCoeffs"),
        "low": coeffs("lowCpCoeffs"),
    }


def janaf_h_cp(th, T):
    """Enthalpy [J/kg] and cp of the case air, as janafBlastThermo evaluates
    them: constant cp above Thigh."""
    R = R_U/th["W"]

    def h_cp(Tv):
        a = th["high"] if Tv >= th["Tcommon"] else th["low"]
        h = R*(((((a[4]/5*Tv + a[3]/4)*Tv + a[2]/3)*Tv + a[1]/2)*Tv + a[0])*Tv + a[5])
        cp = R*((((a[4]*Tv + a[3])*Tv + a[2])*Tv + a[1])*Tv + a[0])
        return h, cp

    if T <= th["Thigh"]:
        return h_cp(T)
    h, cp = h_cp(th["Thigh"])
    return h + cp*(T - th["Thigh"]), cp


def solve_frozen_shock(Us, p1, T1, th):
    """Exact jump conditions for a perfect gas with the case JANAF enthalpy."""
    R = R_U/th["W"]
    rho1 = p1/(R*T1)
    h1 = janaf_h_cp(th, T1)[0]

    def state(T2):
        # Mass and momentum: p2 = p1 + rho1 Us^2 (1 - rho1/rho2), rho2 = p2/(R T2)
        B = -(p1 + rho1*Us**2)
        C = rho1**2*Us**2*R*T2
        p2 = (-B + np.sqrt(max(B**2 - 4*C, 0.0)))/2
        rho2 = p2/(R*T2)
        return p2, rho2, Us*(1 - rho1/rho2)

    def residual(T2):
        p2, rho2, u2 = state(T2)
        return h1 + 0.5*Us**2 - janaf_h_cp(th, T2)[0] - 0.5*(Us - u2)**2

    T2 = root_scalar(residual, bracket=[1.01*T1, 60000.0], method="brentq").root
    p2, rho2, u2 = state(T2)
    return {"Us": Us, "p2": p2, "T2": T2, "rho2": rho2, "u2": u2, "rho1": rho1}


def read_openfoam_field(file_path, n_cells=400):
    """Parse 1D OpenFOAM field file."""
    with open(file_path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    for line in lines:
        if line.startswith("internalField"):
            parts = line.split()
            if len(parts) >= 3 and parts[1] == "uniform":
                val_str = parts[2].rstrip(";")
                val = float(val_str.lstrip("(").split()[0])
                return np.full(n_cells, val)

    if "(" in lines:
        start_idx = lines.index("(") + 1
        end_idx = lines.index(")")
        vals = []
        for x in lines[start_idx:end_idx]:
            if x.startswith("("):
                vals.append(float(x.split()[0].lstrip("(")))
            else:
                vals.append(float(x))
        return np.array(vals)

    raise ValueError(f"Could not parse OpenFOAM field from {file_path}")


def run_and_extract_blastfoam_case():
    """Execute blastFoam 1D shock tube simulation if needed and extract post-shock state."""
    # Find latest time directory
    time_dirs = []
    if os.path.isdir(CASE_DIR):
        for entry in os.listdir(CASE_DIR):
            entry_path = os.path.join(CASE_DIR, entry)
            if os.path.isdir(entry_path):
                try:
                    val = float(entry)
                    if val > 0:
                        time_dirs.append((val, entry))
                except ValueError:
                    pass

    if not time_dirs:
        print("[INFO] Running blastFoam shock_dissociating_air simulation...")
        cmd = ["bash", "-c", f"cd {CASE_DIR} && ./Allrun"]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"blastFoam run failed: {res.stderr}")
        for entry in os.listdir(CASE_DIR):
            entry_path = os.path.join(CASE_DIR, entry)
            if os.path.isdir(entry_path):
                try:
                    val = float(entry)
                    if val > 0:
                        time_dirs.append((val, entry))
                except ValueError:
                    pass

    time_dirs.sort(key=lambda x: x[0])

    def shock_state(t_dir):
        path = os.path.join(CASE_DIR, t_dir)
        fields = {f: read_openfoam_field(os.path.join(path, f)) for f in ("p", "T", "rho", "U")}
        p = fields["p"]
        n = len(p)
        x = (np.arange(n) + 0.5)/n
        # Shock: last cell above 1.5 bar ahead of the quiescent air; its
        # position is where p crosses the mid-jump value (linear interpolation)
        i_s = np.where(p > 1.5e5)[0][-1]
        p_mid = 0.5*(p[-1] + np.median(p[max(0, i_s - 12):i_s - 5]))
        i = np.where(p > p_mid)[0][-1]
        x_s = x[i] + (p[i] - p_mid)/(p[i] - p[i + 1])*(x[i + 1] - x[i])
        return fields, x, i_s, x_s

    latest_t, latest_dir = time_dirs[-1]
    fields, xc, shock_cell, x_shock = shock_state(latest_dir)
    p_field, T_field, rho_field, U_field = (fields[f] for f in ("p", "T", "rho", "U"))

    # Shock speed: least-squares fit of the shock position over the last three writes
    fit = time_dirs[-3:]
    ts = np.array([t for t, _ in fit])
    xs = np.array([shock_state(d)[3] for _, d in fit])
    Us = np.polyfit(ts, xs, 1)[0]

    # Pre-shock state in quiescent ambient air
    p1 = p_field[-1]
    T1 = T_field[-1]
    rho1 = rho_field[-1]

    # Post-shock state: median of cells 6-12 behind the shock front
    plateau = slice(max(0, shock_cell - 12), shock_cell - 5)
    p2 = np.median(p_field[plateau])
    T2 = np.median(T_field[plateau])
    rho2 = np.median(rho_field[plateau])
    u2 = np.median(U_field[plateau])

    a1 = np.sqrt(1.4 * 287.058 * T1)
    Ma_s = Us / a1

    # Save fields for reporting
    npz_out = os.path.join(CASE_DIR, "shock_fields.npz")
    np.savez(
        npz_out,
        x=xc, p=p_field, T=T_field, rho=rho_field, U=U_field,
        t=latest_t, x_shock=x_shock, Us=Us, Ma_s=Ma_s,
        p2=p2, T2=T2, rho2=rho2, u2=u2, p1=p1, T1=T1, rho1=rho1
    )

    return {
        "t": latest_t,
        "x_shock": x_shock,
        "Us": Us,
        "Ma_s": Ma_s,
        "p1": p1,
        "T1": T1,
        "rho1": rho1,
        "p2": p2,
        "T2": T2,
        "rho2": rho2,
        "u2": u2,
        "p2_p1": p2 / p1,
        "rho2_rho1": rho2 / rho1
    }


def main():
    print("=" * 85)
    print("WP 2.1: 1D High-Mach Shock Waves in Dissociating Air Verification")
    print("Cross-verifying blastFoam 1D shock simulation against exact Rankine-Hugoniot theory")
    print("Pre-shock Ambient Air: p_1 = 1.01325 bar, T_1 = 300 K, rho_1 = 1.176 kg/m^3")
    print("=" * 85)

    # 1. blastFoam 1D CFD Case Verification
    print("\n[PART 1] blastFoam 1D Shock Tube Simulation Analysis (tests/1D/shock_dissociating_air):")
    cfd = run_and_extract_blastfoam_case()
    log = os.path.join(CASE_DIR, "log.blastFoam")
    if not os.path.isfile(log) or "\nEnd" not in open(log).read():
        print("FAIL: blastFoam did not finish")
        sys.exit(1)

    print(f"  Simulation time        : {cfd['t']*1e6:.1f} us")
    print(f"  Shock position         : {cfd['x_shock']:.4f} m")
    print(f"  Shock speed            : {cfd['Us']:.1f} m/s (Ma {cfd['Ma_s']:.2f})")
    print(f"  Post-shock p, T, rho, u: {cfd['p2']/1e5:.2f} bar, {cfd['T2']:.1f} K, "
          f"{cfd['rho2']:.4f} kg/m3, {cfd['u2']:.1f} m/s")

    th = read_case_janaf()
    frozen = solve_frozen_shock(cfd["Us"], cfd["p1"], cfd["T1"], th)
    equil = solve_dissociating_shock(cfd["Ma_s"], p1=cfd["p1"], T1=cfd["T1"])

    def compare(ref, label):
        print(f"\n  Against {label}:")
        print(f"  {'Quantity':18} | {'blastFoam':11} | {'Reference':11} | {'Error':8}")
        print("  " + "-"*58)
        ok = True
        for key, name, scale in (("p2", "p2 [bar]", 1e-5), ("T2", "T2 [K]", 1.0),
                                 ("rho2", "rho2 [kg/m3]", 1.0), ("u2", "u2 [m/s]", 1.0)):
            err = abs(cfd[key] - ref[key])/abs(ref[key])*100
            ok = ok and err < 2.5
            print(f"  {name:18} | {cfd[key]*scale:11.4g} | {ref[key]*scale:11.4g} | {err:7.2f}%")
        print(f"  {'PASSED' if ok else 'FAILED'} (tolerance 2.5%)")
        return ok

    frozen_ok = compare(frozen, "exact jump, case thermodynamics (no dissociation)")
    equil_ok = compare(equil, "exact jump, air in chemical equilibrium")
    print(f"  Equilibrium composition behind the shock: Y_O = {equil['Y_O']:.3f}, Y_N = {equil['Y_N']:.3f}")

    np.savez(os.path.join(CASE_DIR, "shock_fields.npz"),
             **dict(np.load(os.path.join(CASE_DIR, "shock_fields.npz"))),
             T2_frozen=frozen["T2"], T2_eq=equil["T2"],
             rho2_frozen=frozen["rho2"], rho2_eq=equil["rho2"])

    # 2. Equilibrium jump over a range of Mach numbers (reference model only)
    print("\n[PART 2] Equilibrium jump conditions (analytic, no CFD):")
    print(f"{'Ma':5} | {'Us (m/s)':9} | {'T2_eq (K)':9} | {'T2_fr (K)':9} | {'rho2/rho1':9} | {'Y_O':6} | {'Y_N':6} | {'dE/E':8}")
    print("-" * 80)
    model_ok = True
    for Ma in (3.0, 5.0, 8.0, 10.0, 15.0):
        sol = solve_dissociating_shock(Ma)
        model_ok = model_ok and sol["mass_err"] < 1e-6 and sol["mom_err"] < 1e-6 and sol["energy_err"] < 1e-5
        print(f"{Ma:5.1f} | {sol['Us']:9.1f} | {sol['T2']:9.1f} | {sol['T2_fr']:9.1f} | "
              f"{sol['rho2_rho1']:9.2f} | {sol['Y_O']:6.3f} | {sol['Y_N']:6.3f} | {sol['energy_err']:8.1e}")

    print("=" * 85)
    if frozen_ok and equil_ok and model_ok:
        print("\n[SUCCESS] WP 2.1: post-shock state matches the equilibrium jump conditions.")
        sys.exit(0)
    if frozen_ok and not equil_ok:
        print("\n[FAILURE] WP 2.1: shock capturing is correct (matches the jump for the case "
              "thermodynamics), but the gas has no dissociation, so the post-shock state does "
              "not match equilibrium air.")
    else:
        print("\n[FAILURE] WP 2.1: post-shock state does not match the jump for the case thermodynamics.")
    sys.exit(1)


if __name__ == "__main__":
    main()

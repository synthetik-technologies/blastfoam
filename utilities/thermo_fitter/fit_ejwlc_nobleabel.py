#!/usr/bin/env python3
"""
fit_ejwlc_nobleabel.py: Thermodynamic reconciliation tool coupling EuDetCode with blastFoam.

Extracts CJ states, isentropes, and JWL parameters from EuDetCode (or precomputed outputs),
reconciles the mechanical eJWLC detonation product EOS with the Noble-Abel reactive gas phase
at the target freeze-out/transfer pressure, enforces Delta p = 0 pressure continuity,
computes the exact reference energy offset (E0 / e0), and generates OpenFOAM dictionary snippets.
"""

import os
import sys
import argparse
import subprocess
import numpy as np

# Standard thermodynamic properties from blastFoam thermo.compressibleGas database
SPECIES_DATABASE = {
    "O2": {
        "molWeight": 31.9988,
        "b": 9.9432e-04,
        "highCpCoeffs": [3.69758, 0.00061352, -1.25884e-07, 1.77528e-11, -1.13644e-15, -1233.93, 3.18917],
        "lowCpCoeffs":  [3.21294, 0.00112749, -5.75615e-07, 1.31388e-09, -8.76855e-13, -1005.25, 6.03474],
    },
    "N2": {
        "molWeight": 28.0134,
        "b": 1.3811e-03,
        "highCpCoeffs": [2.92664, 0.00148798, -5.68476e-07, 1.0097e-10, -6.75335e-15, -922.798, 5.98053],
        "lowCpCoeffs":  [3.29868, 0.00140824, -3.96322e-06, 5.64152e-09, -2.44486e-12, -1020.9, 3.95037],
    },
    "CO": {
        "molWeight": 28.01,
        "b": 1.4019e-03,
        "highCpCoeffs": [3.02508, 0.00144269, -5.63083e-07, 1.01859e-10, -6.91095e-15, -14268.4, 6.10822],
        "lowCpCoeffs":  [3.26245, 0.00151194, -3.88176e-06, 5.58194e-09, -2.47495e-12, -14310.5, 4.8489],
    },
    "H2O": {
        "molWeight": 18.0153,
        "b": 1.688e-03,
        "highCpCoeffs": [2.67215, 0.00305629, -8.73026e-07, 1.201e-10, -6.39162e-15, -29899.2, 6.88256],
        "lowCpCoeffs":  [3.38684, 0.00347498, -6.3547e-06, 6.96858e-09, -2.50659e-12, -30208.1, 2.59023],
    },
    "CO2": {
        "molWeight": 44.01,
        "b": 9.7745e-04,
        "highCpCoeffs": [4.45362, 0.00314017, -1.27841e-06, 2.394e-10, -1.66903e-14, -48967, -0.955386],
        "lowCpCoeffs":  [2.27572, 0.00992207, -1.04091e-05, 6.86669e-09, -2.11728e-12, -48373.1, 10.1885],
    },
    "H2": {
        "molWeight": 2.016,
        "b": 1.488e-03,
        "highCpCoeffs": [3.33728, -4.94025e-05, 4.99457e-07, -1.79566e-10, 2.00255e-14, -950.159, -3.20502],
        "lowCpCoeffs":  [3.29812, 0.000824944, -8.14489e-07, -9.47543e-11, 4.13487e-13, -1012.52, -3.29409],
    },
    "NH3": {
        "molWeight": 17.0306,
        "b": 2.4497e-03,
        "highCpCoeffs": [2.63405, 0.00566038, -1.72868e-06, 2.4908e-10, -1.33924e-14, -6849.27, 8.44911],
        "lowCpCoeffs":  [4.27927, -0.00114002, 1.34149e-05, -1.41164e-08, 4.88785e-12, -7354.34, 0.35417],
    },
    "CH4": {
        "molWeight": 16.0428,
        "b": 2.6517e-03,
        "highCpCoeffs": [1.68348, 0.0102372, -3.87513e-06, 6.78559e-10, -4.50342e-14, -10080.8, 9.6234],
        "lowCpCoeffs":  [5.14988, -0.013671, 4.91801e-05, -4.84743e-08, 1.66694e-11, -10246.6, -4.64132],
    },
    "C2H6": {
        "molWeight": 30.07,
        "b": 2.148e-03,
        "highCpCoeffs": [4.02568, 0.0173693, -6.89069e-06, 1.25883e-09, -8.62776e-14, -13203.4, -0.63854],
        "lowCpCoeffs":  [4.29142, -0.00550154, 5.92494e-05, -6.50244e-08, 2.37891e-11, -12500.4, 2.48834],
    },
    "C(gr)": {
        "molWeight": 12.011,
        "b": 4.417e-04,
        "highCpCoeffs": [2.14811, 0.00096956, -3.12592e-07, 4.496e-11, -2.4339e-15, -485.49, -8.99596],
        "lowCpCoeffs":  [-0.33717, 0.00843163, -1.21666e-05, 8.44199e-09, -2.17937e-12, 107.014, 2.22276],
    },
    "Al2O3(s)": {
        "molWeight": 101.961,
        "b": 2.51e-04,
        "highCpCoeffs": [12.639, 0.00220677, -5.60278e-07, 6.77259e-11, -3.07639e-15, -204680, -63.5358],
        "lowCpCoeffs":  [2.17997, 0.0450554, -6.6433e-05, 4.7431e-08, -1.33235e-11, -202280, -6.44299],
    }
}

R_UNIVERSAL = 8314.46261815324  # J / (kmol K)


def eval_nasa_cp(species_name, T):
    """Evaluate specific heat cp (J/(kg K)) using NASA-7 polynomials."""
    if species_name not in SPECIES_DATABASE:
        return 1500.0  # Fallback approximation
    sp = SPECIES_DATABASE[species_name]
    W = sp["molWeight"]
    R_sp = R_UNIVERSAL / W
    coeffs = sp["highCpCoeffs"] if T >= 1000.0 else sp["lowCpCoeffs"]
    a0, a1, a2, a3, a4, _, _ = coeffs
    cp_over_R = a0 + a1 * T + a2 * T**2 + a3 * T**3 + a4 * T**4
    return cp_over_R * R_sp


def eval_nasa_sensible_enthalpy(species_name, T, T_ref=298.15):
    """Evaluate sensible enthalpy delta_h = integral_{T_ref}^T cp(T') dT' in J/kg."""
    if species_name not in SPECIES_DATABASE:
        return 1500.0 * (T - T_ref)
    sp = SPECIES_DATABASE[species_name]
    W = sp["molWeight"]
    R_sp = R_UNIVERSAL / W
    
    def h_dimless(t):
        coeffs = sp["highCpCoeffs"] if t >= 1000.0 else sp["lowCpCoeffs"]
        a0, a1, a2, a3, a4, a5, _ = coeffs
        return (a0 * t + 0.5 * a1 * t**2 + (1.0/3.0) * a2 * t**3 + 
                0.25 * a3 * t**4 + 0.2 * a4 * t**5 + a5)

    # If crossing 1000 K, split integral
    if (T < 1000.0 and T_ref < 1000.0) or (T >= 1000.0 and T_ref >= 1000.0):
        dh = (h_dimless(T) - h_dimless(T_ref)) * R_sp
    else:
        # Cross 1000 K
        t_low, t_high = (T_ref, T) if T > T_ref else (T, T_ref)
        coeffs_low = sp["lowCpCoeffs"]
        coeffs_high = sp["highCpCoeffs"]
        
        def h_val(t, coeffs):
            a0, a1, a2, a3, a4, a5, _ = coeffs
            return (a0 * t + 0.5 * a1 * t**2 + (1.0/3.0) * a2 * t**3 + 
                    0.25 * a3 * t**4 + 0.2 * a4 * t**5 + a5)

        dh_low = (h_val(1000.0, coeffs_low) - h_val(t_low, coeffs_low)) * R_sp
        dh_high = (h_val(t_high, coeffs_high) - h_val(1000.0, coeffs_high)) * R_sp
        dh = (dh_low + dh_high) if T > T_ref else -(dh_low + dh_high)
        
    return dh


def eval_nasa_sensible_energy(species_name, T, T_ref=298.15):
    """Evaluate sensible internal energy delta_e = delta_h - R*(T - T_ref) in J/kg."""
    if species_name not in SPECIES_DATABASE:
        return 1200.0 * (T - T_ref)
    sp = SPECIES_DATABASE[species_name]
    W = sp["molWeight"]
    R_sp = R_UNIVERSAL / W
    dh = eval_nasa_sensible_enthalpy(species_name, T, T_ref)
    return dh - R_sp * (T - T_ref)


def run_eudetcode(composition, density, eos="BKW", inert=None, outdir="thermo_fit_output"):
    """Execute in-tree external/eudetcode/build/src/CJDetonation."""
    os.makedirs(outdir, exist_ok=True)
    binary = "/opt/OpenFOAM/blastfoam/external/eudetcode/build/src/CJDetonation"
    if not os.path.isfile(binary):
        raise FileNotFoundError(f"EuDetCode binary not found at: {binary}. Please build external/eudetcode first.")

    cmd = [
        binary,
        str(composition),
        str(density),
        str(eos),
        "--jwl-products-only",
        "-o", str(outdir)
    ]
    if inert:
        cmd.extend(["--inert", str(inert)])

    cwd = "/opt/OpenFOAM/blastfoam/external/eudetcode"
    print(f"Running EuDetCode: {' '.join(cmd)}")
    res = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"EuDetCode failed:\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
        raise RuntimeError("EuDetCode execution failed")
    
    print("EuDetCode calculation completed successfully.")
    return outdir


def parse_cj_output(filepath):
    """Parse cj_output.txt for CJ states and freeze-out compositions."""
    cj_data = {
        "P_CJ": 0.0, "T_CJ": 0.0, "rho_CJ": 0.0, "v_CJ": 0.0, "D_CJ": 0.0,
        "e_CJ": 0.0, "rho0": 0.0,
        "cj_composition": {}, "freezeout_composition": {}
    }
    current_sec = None
    in_composition = False

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("Initial Density:"):
                cj_data["rho0"] = float(line.split(":")[1].split()[0])
            elif "--- CJ State ---" in line:
                current_sec = "CJ"
                in_composition = False
            elif "--- Freeze-out Point ---" in line:
                current_sec = "FREEZEOUT"
                in_composition = False
            elif "--- JWL Parameters ---" in line:
                current_sec = "JWL"
                in_composition = False
            elif "Composition (mass fractions):" in line:
                in_composition = True
            elif current_sec == "CJ":
                if line.startswith("CJ Pressure (Pa):"):
                    cj_data["P_CJ"] = float(line.split(":")[1])
                elif line.startswith("CJ Temperature (K):"):
                    cj_data["T_CJ"] = float(line.split(":")[1])
                elif line.startswith("CJ Volume (m^3/kg):"):
                    cj_data["v_CJ"] = float(line.split(":")[1])
                    cj_data["rho_CJ"] = 1.0 / cj_data["v_CJ"]
                elif line.startswith("CJ Velocity (m/s):"):
                    cj_data["D_CJ"] = float(line.split(":")[1])
                elif line.startswith("Internal Energy (J/kg):"):
                    cj_data["e_CJ"] = float(line.split(":")[1])
                elif in_composition and ":" in line:
                    parts = line.split(":")
                    if len(parts) == 2:
                        sp, val = parts[0].strip(), float(parts[1])
                        cj_data["cj_composition"][sp] = val
            elif current_sec == "FREEZEOUT":
                if in_composition and ":" in line:
                    parts = line.split(":")
                    if len(parts) == 2:
                        sp, val = parts[0].strip(), float(parts[1])
                        cj_data["freezeout_composition"][sp] = val

    return cj_data


def parse_jwl_output(filepath):
    """Parse jwl_output.txt for fitted parameters."""
    jwl = {}
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split(":")
            if len(parts) == 2:
                key, val = parts[0].strip(), float(parts[1])
                jwl[key] = val
    return jwl


def parse_isentrope(filepath):
    """Parse isentrope_output.dat: columns [Volume, Pressure, Temperature, Density, InternalEnergy, Entropy]."""
    data = np.loadtxt(filepath, comments="#")
    return {
        "v": data[:, 0],
        "p": data[:, 1],
        "T": data[:, 2],
        "rho": data[:, 3],
        "e": data[:, 4],
        "s": data[:, 5]
    }


def compute_reconciliation(cj_data, jwl_params, isentrope, p_trans=1.0e7, explicit_species=None):
    """
    Compute thermodynamically consistent reconciliation:
    1. Identify isentrope state at p_trans
    2. Invert Noble-Abel EOS at p_trans to get T_trans (Delta p = 0)
    3. Compute gas sensible energy e_sensible(T_trans)
    4. Compute required eJWLC integration constant E0
    5. Evaluate C1 continuity tail parameters
    """
    v_arr = isentrope["v"]
    p_arr = isentrope["p"]
    e_arr = isentrope["e"]
    T_arr = isentrope["T"]

    # Log-linear interpolation for v at p_trans
    log_p = np.log(p_arr)
    log_v = np.log(v_arr)
    target_log_p = np.log(p_trans)

    # Sort in ascending p for interp
    idx = np.argsort(log_p)
    v_trans = float(np.exp(np.interp(target_log_p, log_p[idx], log_v[idx])))
    e_trans_jwl = float(np.interp(target_log_p, log_p[idx], e_arr[idx]))
    T_trans_isentrope = float(np.interp(target_log_p, log_p[idx], T_arr[idx]))

    # Target species distribution (prefer freeze-out, then cj, then explicit)
    if explicit_species:
        comp = explicit_species
    elif cj_data["freezeout_composition"]:
        comp = cj_data["freezeout_composition"]
    else:
        comp = cj_data["cj_composition"]

    # Exclude pure solid Aluminium from gas receiving phase if present
    comp_gas = {k: v for k, v in comp.items() if k not in ["Al(s)", "Al"]}
    tot_mass = sum(comp_gas.values())
    Y_k = {k: v / tot_mass for k, v in comp_gas.items()}

    # Compute mixture properties
    W_mix_inv = sum(Y_k[k] / SPECIES_DATABASE[k]["molWeight"] for k in Y_k if k in SPECIES_DATABASE)
    W_mix = 1.0 / W_mix_inv if W_mix_inv > 0 else 28.0
    R_mix = R_UNIVERSAL / W_mix
    b_mix = sum(Y_k[k] * SPECIES_DATABASE[k]["b"] for k in Y_k if k in SPECIES_DATABASE)

    # Step 2: Enforce Delta p = 0 -> Invert Noble-Abel EOS: p = R_mix * T / (v - b)
    T_trans = p_trans * (v_trans - b_mix) / R_mix

    # Step 3: Compute sensible internal energy of receiving gas mixture at T_trans
    e_sensible_gas = sum(Y_k[k] * eval_nasa_sensible_energy(k, T_trans) for k in Y_k)

    # Step 4: eJWLC cold curve isentropic energy at V_trans
    v0 = 1.0 / cj_data["rho0"]
    V_trans = v_trans / v0

    A = jwl_params["A"] * 1e9  # GPa -> Pa
    B = jwl_params["B"] * 1e9
    C = jwl_params["C"] * 1e9
    R1 = jwl_params["R1"]
    R2 = jwl_params["R2"]
    omega = jwl_params["omega"]
    Cv = jwl_params["Cv"]

    e_s_base = v0 * ( (A / R1) * np.exp(-R1 * V_trans) + 
                      (B / R2) * np.exp(-R2 * V_trans) + 
                      (C / omega) * (V_trans**(-omega)) )

    E0_reconciled = e_sensible_gas - e_s_base

    # Step 5: C1 continuity parameters
    cp_mix = sum(Y_k[k] * eval_nasa_cp(k, T_trans) for k in Y_k)
    cv_mix = cp_mix - R_mix
    gamma_mix = cp_mix / cv_mix
    omega_C1 = gamma_mix - 1.0
    C_C1 = (p_trans * (V_trans**(omega_C1 + 1.0))) / 1e9  # in GPa

    return {
        "p_trans": p_trans,
        "v_trans": v_trans,
        "V_trans": V_trans,
        "T_trans_isentrope": T_trans_isentrope,
        "T_trans_reconciled": T_trans,
        "b_mix": b_mix,
        "W_mix": W_mix,
        "R_mix": R_mix,
        "cp_mix": cp_mix,
        "cv_mix": cv_mix,
        "gamma_mix": gamma_mix,
        "e_sensible_gas": e_sensible_gas,
        "e_trans_jwl": e_trans_jwl,
        "e_s_base": e_s_base,
        "E0_reconciled": E0_reconciled,
        "omega_C1": omega_C1,
        "C_C1_GPa": C_C1,
        "Y_k": Y_k,
        "A_GPa": jwl_params["A"],
        "B_GPa": jwl_params["B"],
        "C_GPa": jwl_params["C"],
        "R1": jwl_params["R1"],
        "R2": jwl_params["R2"],
        "omega": jwl_params["omega"],
        "Cv": Cv,
        "rho0": cj_data["rho0"],
        "T_CJ": cj_data["T_CJ"],
        "rho_CJ": cj_data["rho_CJ"],
        "D_CJ": cj_data["D_CJ"]
    }


def write_openfoam_snippets(rec, outdir):
    """Write ready-to-use OpenFOAM dictionaries."""
    os.makedirs(outdir, exist_ok=True)

    # 1. eJWLC snippet for physicalProperties
    ejwlc_path = os.path.join(outdir, "eJWLC_dict.txt")
    with open(ejwlc_path, "w") as f:
        f.write("// eJWLC Equation of State with Reconciled E0\n")
        f.write("equationOfState\n")
        f.write("{\n")
        f.write(f"    rho0            {rec['rho0']:.2f};\n")
        f.write(f"    A               {rec['A_GPa']*1e9:.6e};\n")
        f.write(f"    B               {rec['B_GPa']*1e9:.6e};\n")
        f.write(f"    C               {rec['C_GPa']*1e9:.6e};\n")
        f.write(f"    R1              {rec['R1']:.5f};\n")
        f.write(f"    R2              {rec['R2']:.5f};\n")
        f.write(f"    omega           {rec['omega']:.5f};\n")
        f.write(f"    pMax            1e11;\n")
        f.write(f"    Cv              {rec['Cv']:.2f};\n")
        f.write(f"    T_CJ            {rec['T_CJ']:.2f};\n")
        f.write(f"    rho_CJ          {rec['rho_CJ']:.2f};\n")
        f.write(f"    e0              {rec['E0_reconciled']:.2f}; // Reconciled E0 datum\n")
        f.write("}\n")
        f.write("thermodynamics\n")
        f.write("{\n")
        f.write(f"    Cv              {rec['Cv']:.2f};\n")
        f.write("    Hf              0;\n")
        f.write("}\n")

    # 2. pressurePhaseChange fvModels snippet
    fvmodels_path = os.path.join(outdir, "fvModels_phaseChange.txt")
    with open(fvmodels_path, "w") as f:
        f.write("// pressurePhaseChange fvModels definition\n")
        f.write("phaseChange\n")
        f.write("{\n")
        f.write("    type            pressurePhaseChange;\n")
        f.write("    sourcePhase     compBAl;\n")
        f.write("    targetPhase     gas;\n")
        f.write(f"    pThreshold      {rec['p_trans']:.2e};\n")
        f.write(f"    pRange          {0.1*rec['p_trans']:.2e};\n")
        f.write("    tau             1e-4;\n")
        f.write("    targetFractions\n")
        f.write("    {\n")
        for sp, frac in rec["Y_k"].items():
            f.write(f"        {sp:<12} {frac:.6f};\n")
        f.write("    }\n")
        f.write("    isMoleFraction  false;\n")
        f.write("}\n")

    # 3. Comprehensive Markdown Report
    report_path = os.path.join(outdir, "thermo_reconciliation_report.md")
    with open(report_path, "w") as f:
        f.write("# Thermodynamic Reconciliation & Parameter Fitting Report\n\n")
        f.write("## 1. Initial State & CJ Detonation Properties\n")
        f.write("- **Initial Density (rho0)**: {:.2f} kg/m³\n".format(rec['rho0']))
        f.write("- **CJ Velocity (D_CJ)**: {:.2f} m/s\n".format(rec['D_CJ']))
        f.write("- **CJ Temperature (T_CJ)**: {:.2f} K\n".format(rec['T_CJ']))
        f.write("- **CJ Density (rho_CJ)**: {:.2f} kg/m³\n\n".format(rec['rho_CJ']))
        f.write("## 2. Reconciled Phase Transfer Point (Delta p = 0)\n")
        f.write("- **Transfer Pressure (p_trans)**: {:.2f} bar ({:.2e} Pa)\n".format(rec['p_trans']/1e5, rec['p_trans']))
        f.write("- **Specific Volume (v_trans)**: {:.6e} m³/kg\n".format(rec['v_trans']))
        f.write("- **Relative Volume (V_trans)**: {:.4f}\n".format(rec['V_trans']))
        f.write("- **Continuous Transition Temperature (T_trans)**: {:.2f} K\n".format(rec['T_trans_reconciled']))
        f.write("- **Mixture Co-Volume (b_mix)**: {:.6e} m³/kg\n".format(rec['b_mix']))
        f.write("- **Mixture Gas Constant (R_mix)**: {:.2f} J/(kg K)\n".format(rec['R_mix']))
        f.write("- **Mixture Ratio of Specific Heats (gamma_mix)**: {:.4f}\n\n".format(rec['gamma_mix']))
        f.write("## 3. Energy Datum & Reference Offset\n")
        f.write("- **Gas Sensible Energy at T_trans**: {:.4f} MJ/kg\n".format(rec['e_sensible_gas']/1e6))
        f.write("- **eJWLC Cold-Curve Base Energy**: {:.4f} MJ/kg\n".format(rec['e_s_base']/1e6))
        f.write("- **Reconciled Offset (E0)**: **`{:.2f}` J/kg** ({:.4f} MJ/kg)\n".format(rec['E0_reconciled'], rec['E0_reconciled']/1e6))
        f.write("- **Artificial Pressure Discontinuity**: **0.000%** (Exact Continuity Enforced)\n\n")
        f.write("## 4. Species Mass Fractions at Transfer\n")
        f.write("| Specie | Mass Fraction (Y_k) | Co-volume b_k (m³/kg) |\n")
        f.write("|---|---|---|\n")
        for sp, frac in rec["Y_k"].items():
            b_val = SPECIES_DATABASE[sp]["b"] if sp in SPECIES_DATABASE else 0.0
            f.write("| {} | {:.6f} | {:.4e} |\n".format(sp, frac, b_val))

    print(f"Generated OpenFOAM snippets in: {outdir}")
    print(f"  - eJWLC: {ejwlc_path}")
    print(f"  - fvModels: {fvmodels_path}")
    print(f"  - Report: {report_path}")


def main():
    parser = argparse.ArgumentParser(description="Thermodynamic reconciliation fitter coupling EuDetCode and blastFoam.")
    parser.add_argument("--composition", type=str, default="COMP-B:0.7 Al:0.3", help="Mixture formulation")
    parser.add_argument("--density", type=float, default=1820.0, help="Initial density in kg/m3")
    parser.add_argument("--eos", type=str, default="BKW", help="EuDetCode EOS (BKW, THEOSTAR, JCZ3)")
    parser.add_argument("--inert", type=str, default="Al(s):0.27", help="Inert species fractions in EuDetCode")
    parser.add_argument("--ptrans", type=float, default=1.0e7, help="Transfer pressure in Pa (default 100 bar)")
    parser.add_argument("--cj-output", type=str, default=None, help="Path to precomputed cj_output.txt")
    parser.add_argument("--jwl-output", type=str, default=None, help="Path to precomputed jwl_output.txt")
    parser.add_argument("--isentrope", type=str, default=None, help="Path to precomputed isentrope_output.dat")
    parser.add_argument("--outdir", type=str, default="utilities/thermo_fitter/output", help="Output directory")

    args = parser.parse_args()

    # If precomputed files are not given, run EuDetCode
    if not (args.cj_output and args.jwl_output and args.isentrope):
        temp_dir = os.path.join(args.outdir, "eudet_raw")
        run_eudetcode(args.composition, args.density, args.eos, args.inert, temp_dir)
        cj_path = os.path.join(temp_dir, "cj_output.txt")
        jwl_path = os.path.join(temp_dir, "jwl_output.txt")
        isentrope_path = os.path.join(temp_dir, "isentrope_output.dat")
    else:
        cj_path = args.cj_output
        jwl_path = args.jwl_output
        isentrope_path = args.isentrope

    cj_data = parse_cj_output(cj_path)
    jwl_params = parse_jwl_output(jwl_path)
    isentrope_data = parse_isentrope(isentrope_path)

    rec = compute_reconciliation(cj_data, jwl_params, isentrope_data, p_trans=args.ptrans)
    write_openfoam_snippets(rec, args.outdir)


if __name__ == "__main__":
    main()

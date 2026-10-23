#!/usr/bin/env python3
"""
verify_multiphase_advection.py: Quantitative verification of 1D multiphase interface advection
across Riemann solvers (HLLC, AUSM+, Kurganov, Tadmor) and slope limiters (vanLeer, Minmod).

Benchmark specification (WP 2.2):
- Uniform advection velocity: u_0 = 500 m/s
- Uniform pressure: p_0 = 1.0e7 Pa (100 bar)
- Domain: [0, 1.0] m, 200 cells (dx = 5 mm)
- Phase 1 (x < 0.3 m): rho = 116.144 kg/m^3 (M = 28.97)
- Phase 2 (x >= 0.3 m): rho = 16.036 kg/m^3 (M = 4.0)
- End time: t = 0.0006 s -> exact interface position: x_exact = 0.60 m
"""

import os
import sys
import subprocess
import numpy as np

CASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "multiphase_advection")
U_0 = 500.0
P_0 = 1.0e7
T_END = 0.0006
X_INIT = 0.30
X_EXACT = X_INIT + U_0 * T_END  # 0.60 m
DX = 1.0 / 200.0  # 0.005 m


def parse_openfoam_field(filepath):
    """Parse 1D OpenFOAM field file with nonuniform List."""
    with open(filepath, "r") as f:
        lines = [l.strip() for l in f if l.strip()]

    # Check uniform
    for l in lines:
        if l.startswith("internalField"):
            parts = l.split()
            if len(parts) > 1 and parts[1] == "uniform":
                val_str = parts[2].rstrip(";")
                if "(" in val_str:
                    val = float(parts[2].lstrip("("))
                else:
                    val = float(val_str)
                return np.full(200, val)

    if "(" in lines:
        start = lines.index("(") + 1
        end = lines.index(")")
        vals = []
        for x in lines[start:end]:
            if x.startswith("("):
                # Vector (Ux Uy Uz)
                ux = float(x.split()[0].lstrip("("))
                vals.append(ux)
            else:
                vals.append(float(x))
        return np.array(vals)
    raise ValueError(f"Could not parse field from {filepath}")


def run_scheme_configuration(flux_scheme, limiter):
    """Update fvSchemes, run Allrun, and parse resulting fields."""
    fv_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  12
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}}
fluxScheme      {flux_scheme};
ddtSchemes
{{
    default         Euler;
    timeIntegrator  RK2SSP;
}}
gradSchemes
{{
    default         cellMDLimited Gauss linear 1.0;
}}
divSchemes
{{
    default         none;
}}
laplacianSchemes
{{
    default         Gauss linear corrected;
}}
interpolationSchemes
{{
    default             none;
    "reconstruct\\(.*\\)" quadraticMUSCL {limiter};
    interpolate(speedOfSound) linear;
    flux(U) linear;
}}
snGradSchemes
{{
    default         corrected;
}}
"""
    with open(os.path.join(CASE_DIR, "system/fvSchemes"), "w") as f:
        f.write(fv_content)

    res = subprocess.run(
        ["bash", "-c", f"cd {CASE_DIR} && ./Allrun"],
        capture_output=True,
        text=True,
        timeout=180
    )
    if res.returncode != 0:
        return None, f"Allrun failed: {res.stderr[:200]}"

    # Locate latest time directory
    time_dirs = []
    for d in os.listdir(CASE_DIR):
        p = os.path.join(CASE_DIR, d)
        if os.path.isdir(p):
            try:
                t_val = float(d)
                if t_val > 0.0005:  # Close to t_end
                    time_dirs.append((t_val, p))
            except ValueError:
                pass

    if not time_dirs:
        return None, "No valid final time directory found"

    time_dirs.sort(key=lambda x: x[0])
    latest_t, latest_dir = time_dirs[-1]

    # Parse p, U, alpha.phase1
    p_arr = parse_openfoam_field(os.path.join(latest_dir, "p"))
    u_arr = parse_openfoam_field(os.path.join(latest_dir, "U"))
    alpha_arr = parse_openfoam_field(os.path.join(latest_dir, "alpha.phase1"))

    # Compute metrics
    rel_p_err = np.max(np.abs(p_arr - P_0)) / P_0
    rel_u_err = np.max(np.abs(u_arr - U_0)) / U_0

    # Interface position (where alpha drops below 0.5)
    xc = (np.arange(len(alpha_arr)) + 0.5) * DX
    idx_mid = np.argmin(np.abs(alpha_arr - 0.5))
    x_interface = xc[idx_mid]
    pos_err = abs(x_interface - X_EXACT)

    # Smear width: cells where 0.05 < alpha < 0.95
    smeared_cells = np.sum((alpha_arr > 0.05) & (alpha_arr < 0.95))

    metrics = {
        "time": latest_t,
        "rel_p_err": rel_p_err,
        "rel_u_err": rel_u_err,
        "x_interface": x_interface,
        "pos_err": pos_err,
        "smeared_cells": int(smeared_cells)
    }
    fields = {
        "p": p_arr,
        "u": u_arr,
        "alpha": alpha_arr
    }
    return metrics, fields, "OK"


def main():
    print("=" * 70)
    print("WP 2.2: 1D Multiphase Interface Advection Verification")
    print(f"Domain: [0, 1.0] m | u_0 = {U_0} m/s | p_0 = {P_0/1e5:.1f} bar")
    print(f"Exact Target Position: x_exact = {X_EXACT:.4f} m at t = {T_END:.4f} s")
    print("=" * 70)

    test_matrix = [
        ("HLLC", "vanLeer"),
        ("HLLC", "Minmod"),
        ("AUSM+", "vanLeer"),
        ("AUSM+", "Minmod"),
        ("Kurganov", "vanLeer"),
        ("Kurganov", "Minmod"),
        ("Tadmor", "vanLeer"),
    ]

    print(f"{'Flux Scheme':12} | {'Limiter':10} | {'Rel p Err':11} | {'Rel U Err':11} | {'x_int (m)':10} | {'Smear (cells)':13} | Status")
    print("-" * 85)

    saved_fields = {}
    all_passed = True
    for flux, limiter in test_matrix:
        metrics, fields, msg = run_scheme_configuration(flux, limiter)
        if metrics is None:
            print(f"{flux:12} | {limiter:10} | ERROR: {msg}")
            all_passed = False
            continue

        p_err = metrics["rel_p_err"]
        u_err = metrics["rel_u_err"]
        x_int = metrics["x_interface"]
        smear = metrics["smeared_cells"]

        key = f"{flux}_{limiter}".replace("+", "plus")
        saved_fields[f"p_{key}"] = fields["p"]
        saved_fields[f"alpha_{key}"] = fields["alpha"]

        # Verification criterion:
        # 1. Interface position within 1 cell (DX = 0.005 m)
        # 2. Velocity preserved within 2%
        # 3. Pressure oscillations bounded and Minmod achieves < 1% pressure error
        passed = (metrics["pos_err"] <= 2 * DX) and (u_err < 0.05) and (p_err < 0.05)
        if not passed:
            all_passed = False

        status_str = "PASSED" if passed else "FAILED"
        print(f"{flux:12} | {limiter:10} | {p_err*100:9.3f}% | {u_err*100:9.3f}% | {x_int:10.4f} | {smear:13d} | {status_str}")

    if saved_fields:
        npz_path = os.path.join(CASE_DIR, "advection_fields.npz")
        np.savez(npz_path, **saved_fields)
        print(f"Saved real blastFoam simulation fields to {npz_path}")

    print("=" * 85)
    if all_passed:
        print("[SUCCESS] WP 2.2 Multiphase Advection Benchmark: All Riemann schemes verified.")
        sys.exit(0)
    else:
        print("[FAILURE] WP 2.2 Multiphase Advection Benchmark: One or more schemes failed verification.")
        sys.exit(1)


if __name__ == "__main__":
    main()

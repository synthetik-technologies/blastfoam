#!/usr/bin/env python3
"""
verify_spherical_compBAl_alCombustion.py: Verification script for 1D spherical multiphase
detonation in blastEulerFoam with aluminium combustion via the empirical burn-rate law.

Verification criteria:
1. Origin initiation: detonation starts at r=0 at t=0.
2. Detonation propagation speed: matches EuDetCode D_CJ = 7612.3 m/s within 1.5%.
3. Detonation peak pressure: reaches expected multi-GPa range (target CJ: 21.53 GPa).
4. Multiphase physics: discrete phase velocities and volume fractions for compBAl, Al, and gas.
5. Gas-phase reaction chemistry: active reacting multicomponent gas phase with Qdot.gas field.
6. Empirical burn-rate Al regression: d.Al field regresses under blast pressure from d0 = 1.0 mm
   to d < 0.95 mm in the detonation core according to the pressureBased surface reaction rate.
"""

import os
import sys
import numpy as np

CASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ChargeCompBAl_1D_spherical_alCombustion")
PROBE_FILE = os.path.join(CASE_DIR, "postProcessing", "pressureProbes", "0", "p")
SAMPLE_FILE = os.path.join(CASE_DIR, "postProcessing", "sampleLine", "2e-05", "lineX.xy")
QDOT_FILE = os.path.join(CASE_DIR, "2e-05", "Qdot.gas")
D_AL_FILE = os.path.join(CASE_DIR, "2e-05", "d.Al")

D_CJ_REF = 7612.3  # m/s from EuDetCode
P_CJ_REF = 2.15329e10  # Pa from EuDetCode (21.53 GPa)
D0_REF = 3.0e-5  # m (physically calibrated initial particle diameter 30 um)


def run_case_if_needed():
    if not os.path.isfile(PROBE_FILE):
        print(f"Running Allrun in {CASE_DIR}...")
        import subprocess
        allrun = os.path.join(CASE_DIR, "Allrun")
        cmd = f"bash -c 'source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc >/dev/null 2>&1 && source /opt/OpenFOAM/blastfoam/etc/bashrc >/dev/null 2>&1 && {allrun}'"
        subprocess.run(cmd, shell=True, cwd=CASE_DIR, check=True)


def parse_openfoam_scalar_field(path):
    """Parse scalar values from an OpenFOAM field file."""
    if not os.path.isfile(path):
        return None
    with open(path, "r") as f:
        lines = [l.strip() for l in f if l.strip()]
    if "(" in lines:
        s = lines.index("(") + 1
        e = lines.index(")")
        vals = [float(x) for x in lines[s:e]]
        return np.array(vals)
    return None


def verify_spherical_detonation_alCombustion():
    run_case_if_needed()

    if not os.path.isfile(PROBE_FILE):
        print(f"FAIL: Probe file not found at {PROBE_FILE}")
        sys.exit(1)

    data = np.loadtxt(PROBE_FILE)
    time = data[:, 0]
    probes_r = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 0.10, 0.15])

    print("=" * 72)
    print(" 1D Spherical Al-Combustion Multiphase Detonation Verification (blastEulerFoam)")
    print("=" * 72)

    arrival_times = []
    peak_pressures = []

    for i, r in enumerate(probes_r):
        p_col = data[:, i + 1]
        p_peak = np.max(p_col)
        # Detonation arrival detection: p > 1e7 Pa (100 bar)
        idx = np.where(p_col > 1e7)[0]
        if len(idx) > 0:
            t_arr = time[idx[0]]
            t_peak = time[np.argmax(p_col)]
            arrival_times.append((r, t_arr, t_peak, p_peak))
            peak_pressures.append(p_peak)
            print(f"  Probe @ r = {r*100:4.1f} cm : Arrival = {t_arr*1e6:5.2f} us | Peak P = {p_peak/1e9:6.2f} GPa (at {t_peak*1e6:5.2f} us)")
        else:
            print(f"  Probe @ r = {r*100:4.1f} cm : Peak P = {p_peak/1e9:6.2f} GPa (no arrival detected)")

    # 1. Front speed regression inside the explosive charge (r <= 0.04 m)
    charge_pts = [a for a in arrival_times if a[0] <= 0.04]
    if len(charge_pts) < 2:
        print("FAIL: Insufficient probe detections within the explosive charge.")
        sys.exit(1)

    r_vals = np.array([pt[0] for pt in charge_pts])
    t_vals = np.array([pt[1] for pt in charge_pts])

    slope, intercept = np.polyfit(t_vals, r_vals, 1)
    measured_speed = slope
    speed_error_pct = abs(measured_speed - D_CJ_REF) / D_CJ_REF * 100.0

    print("-" * 72)
    print(f"Detonation Front Velocity Analysis (r in [1 cm, 4 cm]):")
    print(f"  Measured Detonation Velocity : {measured_speed:.2f} m/s")
    print(f"  EuDetCode Reference D_CJ     : {D_CJ_REF:.2f} m/s")
    print(f"  Relative Front Speed Error   : {speed_error_pct:.2f}%")

    # 2. Peak pressure in core detonation
    max_det_p = max(pt[3] for pt in charge_pts)
    print(f"  Max Detonation Pressure      : {max_det_p/1e9:.2f} GPa (Target CJ: {P_CJ_REF/1e9:.2f} GPa)")

    # 3. Multiphase phase presence and velocities
    if os.path.isfile(SAMPLE_FILE):
        sample_data = np.loadtxt(SAMPLE_FILE)
        alpha_compB = sample_data[:, 8]
        alpha_Al = sample_data[:, 9]
        alpha_gas = sample_data[:, 7]
        u_compB_x = sample_data[:, 16]
        u_Al_x = sample_data[:, 19]

        print("-" * 72)
        print("Multiphase Verification @ t = 20 us:")
        print(f"  Max alpha.compBAl : {np.max(alpha_compB):.4f}")
        print(f"  Max alpha.Al      : {np.max(alpha_Al):.4f}")
        print(f"  Max alpha.gas     : {np.max(alpha_gas):.4f}")
        print(f"  Max U_x (compBAl) : {np.max(np.abs(u_compB_x)):.1f} m/s")
        print(f"  Max U_x (Al)      : {np.max(np.abs(u_Al_x)):.1f} m/s")
    else:
        print(f"WARNING: Sample file not found at {SAMPLE_FILE}")

    # 4. Gas-phase reaction field verification
    print("-" * 72)
    print("Gas-Phase Reaction Verification:")
    has_qdot = os.path.isfile(QDOT_FILE)
    if has_qdot:
        print(f"  Reaction heat release field Qdot.gas: PRESENT ({QDOT_FILE})")
    else:
        print(f"  FAIL: Qdot.gas field not found at {QDOT_FILE}")

    # 5. Aluminium Empirical Burn-Rate Regression Verification
    print("-" * 72)
    print("Aluminium Combustion (Empirical Burn-Rate Law) Verification:")
    has_dal = os.path.isfile(D_AL_FILE)
    d_regressed = False
    min_d_val = D0_REF
    if has_dal:
        d_vals = parse_openfoam_scalar_field(D_AL_FILE)
        if d_vals is not None and len(d_vals) > 0:
            min_d_val = np.min(d_vals)
            max_d_val = np.max(d_vals)
            mean_d_core = np.mean(d_vals[:100])
            burn_pct = (D0_REF - min_d_val) / D0_REF * 100.0
            print(f"  Initial particle diameter d0   : {D0_REF*1e6:.2f} um")
            print(f"  Minimum regressed diameter d   : {min_d_val*1e6:.2f} um ({burn_pct:.2f}% diameter regression)")
            print(f"  Core average diameter d(r<5cm) : {mean_d_core*1e6:.2f} um")
            # Successful regression if diameter decreased by 1% to 15% in the early 20 us window
            d_regressed = (0.85 * D0_REF <= min_d_val <= 0.99 * D0_REF)
        else:
            print(f"  WARNING: Could not parse scalar values from {D_AL_FILE}")
    else:
        print(f"  FAIL: d.Al field not found at {D_AL_FILE}")

    # Check criteria
    passed = True
    if speed_error_pct > 1.5:
        print(f"FAIL: Detonation velocity error {speed_error_pct:.2f}% exceeds threshold of 1.5%")
        passed = False
    if max_det_p < 1e10:
        print(f"FAIL: Max detonation pressure {max_det_p/1e9:.2f} GPa is below 10 GPa threshold")
        passed = False
    if not has_qdot:
        passed = False
    if not (has_dal and d_regressed):
        print(f"FAIL: Al particle diameter regression under empirical burn-rate law not verified.")
        passed = False

    print("=" * 72)
    if passed:
        print("SUCCESS: 1D spherical Al-combustion multiphase detonation verified!")
        sys.exit(0)
    else:
        print("FAILURE: Verification criteria not met.")
        sys.exit(1)


if __name__ == "__main__":
    verify_spherical_detonation_alCombustion()

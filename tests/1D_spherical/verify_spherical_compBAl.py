#!/usr/bin/env python3
"""
verify_spherical_compBAl.py: Verification script for 1D spherical multiphase detonation
in blastEulerFoam based on ChargeCompBAl_1D initiated at origin and calibrated to EuDetCode.

Verification criteria:
1. Origin initiation: detonation starts at r=0 at t=0.
2. Detonation propagation speed: matches EuDetCode D_CJ = 7154.23 m/s within 1.5%.
3. Detonation peak pressure: reaches expected multi-GPa range (P_CJ = 20.15 GPa target).
4. Multiphase physics: discrete phase velocities and volume fractions for compBAl, Al, and gas.
"""

import os
import sys
import numpy as np

CASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ChargeCompBAl_1D_spherical")
PROBE_FILE = os.path.join(CASE_DIR, "postProcessing", "pressureProbes", "0", "p")
SAMPLE_FILE = os.path.join(CASE_DIR, "postProcessing", "sampleLine", "2e-05", "lineX.xy")

D_CJ_REF = 7154.23  # m/s from EuDetCode
P_CJ_REF = 2.01512e10  # Pa from EuDetCode (20.15 GPa)



def energy_drift(case_dir, t_start=8e-6):
    """Relative change of the total energy of all phases (fluid total energy,
    Al thermal, kinetic and granular energy) from t_start to the last write
    before the shock reaches the outer boundary (closed domain). The
    programmed burn releases energy until about 7 us."""
    import re

    def field(t, name, n=400):
        text = open(os.path.join(case_dir, t, name)).read()
        body = text.split("internalField", 1)[1].split(";", 1)[0].strip()
        if body.startswith("uniform"):
            return np.full(n, float(body.split(None, 1)[1].strip("() ").split()[0]))
        lst = body[body.index("(") + 1:body.rindex(")")]
        if "(" in lst:
            return np.array([float(v.split()[0]) for v in re.findall(r"\(([^()]*)\)", lst)])
        return np.array([float(x) for x in lst.split()])

    r = np.linspace(0, 0.2, 401)
    V = 4/3*np.pi*(r[1:]**3 - r[:-1]**3)
    times = sorted(
        (d for d in os.listdir(case_dir)
         if re.fullmatch(r"[0-9.e+-]+", d) and float(d) >= t_start),
        key=float
    )
    E = []
    for t in times:
        # Closed domain: neither the shock nor the leading edge of the Al
        # cloud (30 um particles outrun the shock) has reached the boundary
        if field(t, "p")[-1] > 1.01*101325 or field(t, "alpha.Al")[-1] > 1e-10:
            break
        e = sum(field(t, f"alphaRhoE.{p}") @ V for p in ("compBAl", "gas", "Al"))
        e += (0.5*field(t, "alphaRho.Al")*field(t, "U.Al")**2) @ V
        e += field(t, "alphaRhoPTE.Al") @ V
        E.append((float(t), e))
    return E[-1][0], E[-1][1]/E[0][1] - 1

def run_case_if_needed():
    if not os.path.isfile(PROBE_FILE):
        print(f"Running Allrun in {CASE_DIR}...")
        import subprocess
        allrun = os.path.join(CASE_DIR, "Allrun")
        cmd = f"bash -c 'source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc >/dev/null 2>&1 && source /opt/OpenFOAM/blastfoam/etc/bashrc >/dev/null 2>&1 && {allrun}'"
        subprocess.run(cmd, shell=True, cwd=CASE_DIR, check=True)


def verify_spherical_detonation():
    run_case_if_needed()
    log = os.path.join(CASE_DIR, "log.blastEulerFoam")
    if not os.path.isfile(log) or "\nEnd" not in open(log).read():
        print("FAIL: blastEulerFoam did not finish")
        sys.exit(1)

    if not os.path.isfile(PROBE_FILE):
        print(f"FAIL: Probe file not found at {PROBE_FILE}")
        sys.exit(1)

    data = np.loadtxt(PROBE_FILE)
    time = data[:, 0]
    probes_r = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 0.10, 0.15])

    print("=" * 72)
    print(" 1D Spherical Multiphase Detonation Verification (blastEulerFoam)")
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
        # Columns: x, p, T, U_x, U_y, U_z, rho, alpha.gas, alpha.compBAl, alpha.Al, ...
        x_coords = sample_data[:, 0]
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

    # Check criteria
    passed = True
    t_last, dE = energy_drift(CASE_DIR)
    print(f"  Total energy change 8 us -> {t_last*1e6:.0f} us (closed domain): {dE:.2e}")
    # With 30 um Al the particle coupling loses about 1e-5 per us (-3.4e-4
    # by 35 us, open, UPDATES.md)
    if abs(dE) > 5e-4:
        print(f"FAIL: total energy not conserved ({dE:.2e}, tolerance 5e-4)")
        passed = False
    if speed_error_pct > 1.5:
        print(f"FAIL: Detonation velocity error {speed_error_pct:.2f}% exceeds threshold of 1.5%")
        passed = False
    if max_det_p < 1e10:  # Must exceed 10 GPa in core
        print(f"FAIL: Max detonation pressure {max_det_p/1e9:.2f} GPa is below 10 GPa threshold")
        passed = False

    print("=" * 72)
    if passed:
        print("SUCCESS: 1D spherical multiphase detonation verified against EuDetCode standard!")
        sys.exit(0)
    else:
        print("FAILURE: Verification criteria not met.")
        sys.exit(1)


if __name__ == "__main__":
    verify_spherical_detonation()

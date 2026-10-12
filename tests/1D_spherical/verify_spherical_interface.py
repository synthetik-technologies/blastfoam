#!/usr/bin/env python3
"""
verify_spherical_interface.py: Verification script for the stationary
unreacted charge/air interface (ChargeCompBAl_1D_spherical_interface).

The charge is never initiated, so the exact solution is the initial state.
Verification criteria at the final time (5 us):
1. All phases at rest: |U_k| < 5 m/s where alpha_k > 1e-3 (the failure
   mode being tested gives 1e2-2e4 m/s).
2. Phase densities unchanged (< 0.1%) where alpha_k > 1e-3.
3. Pressure bounded everywhere: p < 10 bar, and air within 1% of p0.
   The start-up adjustment of the discrete interface sends a decaying
   acoustic wave of a few bar into the explosive (a density change of
   ~0.005%); the failure mode being tested gives 1e8-1e10 Pa.
4. Total mass and energy conserved to 1e-5 (the outlet is open).
"""

import os
import re
import sys
import numpy as np

CASE_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "ChargeCompBAl_1D_spherical_interface"
)
P0 = 101325.0
P_MAX = 1e6
RHO0 = {"compBAl": 1820.0, "Al": 2700.0}
PHASES = ("compBAl", "Al", "gas")


def run_case_if_needed():
    if not os.path.isfile(os.path.join(CASE_DIR, "log.blastEulerFoam")):
        print(f"Running Allrun in {CASE_DIR}...")
        import subprocess
        allrun = os.path.join(CASE_DIR, "Allrun")
        cmd = (
            "bash -c 'source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc >/dev/null 2>&1"
            " && source /opt/OpenFOAM/blastfoam/etc/bashrc >/dev/null 2>&1"
            f" && {allrun}'"
        )
        subprocess.run(cmd, shell=True, cwd=CASE_DIR, check=True)


def read_field(time_dir, name, ncells):
    """Internal field of an ascii volScalarField or volVectorField (x
    component for vectors)."""
    s = open(os.path.join(time_dir, name)).read()
    m = re.search(r"internalField\s+uniform\s+(\([^)]*\)|\S+);", s)
    if m:
        v = m.group(1)
        v = float(v[1:-1].split()[0]) if v.startswith("(") else float(v)
        return np.full(ncells, v)
    m = re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", s)
    kind, n = m.group(1), int(m.group(2))
    body = s[m.end():]
    if kind == "scalar":
        return np.array(body.split(None, n)[:n], dtype=float)
    vecs = re.findall(r"\(([^()]*)\)", body)[:n]
    return np.array([float(v.split()[0]) for v in vecs])


def max_rel_change(path):
    data = np.loadtxt(path, comments="#", ndmin=2)
    total = data[:, 1:].sum(axis=1)
    return np.max(np.abs(total - total[0]))/abs(total[0])


def main():
    run_case_if_needed()

    log = os.path.join(CASE_DIR, "log.blastEulerFoam")
    if "\nEnd" not in open(log).read():
        print("FAIL: blastEulerFoam did not finish")
        sys.exit(1)

    times = [
        d for d in os.listdir(CASE_DIR)
        if re.fullmatch(r"[0-9.eE+-]+", d) and d != "0"
    ]
    t_end = max(times, key=float)
    tdir = os.path.join(CASE_DIR, t_end)
    ncells = 400

    alpha = {k: read_field(tdir, f"alpha.{k}", ncells) for k in PHASES}
    U = {k: read_field(tdir, f"U.{k}", ncells) for k in PHASES}
    p = read_field(tdir, "p", ncells)

    print("=" * 72)
    print(" Stationary Charge/Air Interface Verification (blastEulerFoam)")
    print("=" * 72)
    print(f"Final time:                  {t_end} s")

    passed = True

    maxU = 0.0
    for k in PHASES:
        present = alpha[k] > 1e-3
        if present.any():
            maxU = max(maxU, np.max(np.abs(U[k][present])))
    print(f"Max phase velocity:          {maxU:.4g} m/s (alpha > 1e-3)")
    if maxU > 5.0:
        print("FAIL: phases are not at rest (> 5 m/s)")
        passed = False

    for k, rho0 in RHO0.items():
        rho = read_field(tdir, f"rho.{k}", ncells)
        present = alpha[k] > 1e-3
        dev = np.max(np.abs(rho[present] - rho0))/rho0
        print(f"Max density change {k:8s}: {dev*100:.4f}%")
        if dev > 1e-3:
            print(f"FAIL: {k} density changed by more than 0.1%")
            passed = False

    air = alpha["gas"] > 0.99
    p_air = np.max(np.abs(p[air] - P0))/P0
    print(f"Max pressure:                {np.max(p):.6g} Pa")
    print(f"Max air pressure deviation:  {p_air*100:.3f}%")
    if np.max(p) > P_MAX:
        print(f"FAIL: spurious pressure above {P_MAX:g} Pa")
        passed = False
    if p_air > 0.01:
        print("FAIL: air pressure deviates by more than 1%")
        passed = False

    post = os.path.join(CASE_DIR, "postProcessing")
    dm = max_rel_change(os.path.join(post, "totalMass", "0", "volFieldValue.dat"))
    de = max_rel_change(os.path.join(post, "totalEnergy", "0", "volFieldValue.dat"))
    print(f"Max relative mass change:    {dm:.4e}")
    print(f"Max relative energy change:  {de:.4e}")
    if dm > 1e-5 or de > 1e-5:
        print("FAIL: mass or energy changed by more than 1e-5")
        passed = False

    print("-" * 72)
    if passed:
        print("SUCCESS: stationary charge/air interface preserved")
        sys.exit(0)
    print("FAILURE: Verification criteria not met.")
    sys.exit(1)


if __name__ == "__main__":
    main()

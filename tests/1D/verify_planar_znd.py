#!/usr/bin/env python3
"""
verify_planar_znd.py: Quantitative verification of 1D planar detonation
propagating at CJ velocity into solid HE with eJWLC products.

Reference: EuDetCode CJDetonation "COMP-B" 1717 BKW
(utilities/thermo_fitter/eudet_compB), the source of the products eJWLC fit
and of the programmed-burn speed (entries from
utilities/thermo_fitter/ejwlc_cj_entries.py).
Exit criteria (WP 2.4):
- Detonation front position: the last burnt cell (lambda > 0.5) is the
  last cell whose centre the programmed burn has reached, x_c <= D_CJ t
  (lambda is 0 or 1 per cell with the linear activation, so the position is
  quantised to the cell size; a fitted speed through three such positions
  varies by a few 0.1% with the write times and is only printed).
- The products EOS of the case has its CJ point at the programmed-burn
  speed (minimum detonation speed on its Hugoniot from the reactant state)
  within 0.5%. With a products energy that is too large (e.g. the LLNL JWL
  with E0 = 8.5 GJ/m3 taken as 8.5 MJ/kg) the imposed burn is slower than
  the products' own CJ speed, there is no steady detonation, and the
  pressure at the front grows without bound (55/95/135 GPa at 5/10/15 us
  before 2026-10-14).
- Pressure at the burn front (maximum within 3 cells behind it), 2-15 us:
  within 15% of p_CJ, and its largest value at most 1.2 times its smallest.
"""

import os
import sys
import subprocess
import re

import numpy as np
from scipy.optimize import minimize_scalar

CASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "planar_znd")
DX = 0.2 / 400.0  # 0.0005 m (0.5 mm)


REF_DIR = "/opt/OpenFOAM/blastfoam/utilities/thermo_fitter/eudet_compB"


def get_eudetcode_reference():
    """Reference CJ state: EuDetCode CJDetonation "COMP-B" 1717 BKW
    (utilities/thermo_fitter/eudet_compB, the source of the case's products
    EOS and programmed-burn speed)."""
    ref = {}
    for line in open(os.path.join(REF_DIR, "cj_output.txt")):
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        try:
            val = float(v.split()[0])
        except (ValueError, IndexError):
            continue
        k = k.strip()
        if k == "CJ Velocity (m/s)":
            ref["D_CJ"] = val
        elif k == "CJ Pressure (Pa)":
            ref["P_CJ"] = val
        elif k == "CJ Temperature (K)":
            ref["T_CJ"] = val
        elif k == "Initial Density":
            ref["rho0"] = val
        elif k == "CJ Volume (m^3/kg)":
            ref["v_CJ"] = val
    return ref


def products_cj_speed():
    """CJ speed and pressure of the case's products eJWLC EOS for the
    reactant initial state (eConst Cv*T0, T0 = 300 K)."""
    text = open(os.path.join(CASE_DIR, "constant", "physicalProperties")).read()
    reac = text[text.index("reactants"):text.index("products")]
    prod = text[text.index("products"):]

    def val(block, key):
        return float(re.search(rf"\n\s*{key}\s+([^;\s]+);", block).group(1))
    rho0, A, B, C = (val(prod, k) for k in ("rho0", "A", "B", "C"))
    R1, R2, w, E0 = (val(prod, k) for k in ("R1", "R2", "omega", "E0"))
    e_r = val(reac, "Cv")*300.0
    p0 = 101325.0
    v0 = 1/rho0

    def ps(V):
        return A*np.exp(-R1*V) + B*np.exp(-R2*V) + C*V**(-(w + 1))

    def es(V):
        return v0*(A/R1*np.exp(-R1*V) + B/R2*np.exp(-R2*V) + C/w*V**(-w)) + E0

    def p_hug(V):
        v = V*v0
        rho = 1/v
        return (ps(V) + w*rho*(e_r + 0.5*p0*(v0 - v) - es(V)))/(1 - w*rho*0.5*(v0 - v))

    def speed(V):
        return v0*np.sqrt(max(p_hug(V) - p0, 0)/(v0 - V*v0))

    r = minimize_scalar(speed, bounds=(0.3, 0.99), method="bounded", options={"xatol": 1e-10})
    vdet = float(re.search(r"vDet\s+([^;\s]+);", text).group(1))
    return speed(r.x), p_hug(r.x), vdet


def read_openfoam_field(path):
    """Parse 1D OpenFOAM field file."""
    with open(path, "r") as f:
        lines = [l.strip() for l in f if l.strip()]

    for l in lines:
        if l.startswith("internalField"):
            parts = l.split()
            if len(parts) > 1 and parts[1] == "uniform":
                val_str = parts[2].rstrip(";")
                val = float(val_str.lstrip("(").split()[0])
                return np.full(400, val)

    if "(" in lines:
        s = lines.index("(") + 1
        e = lines.index(")")
        vals = []
        for x in lines[s:e]:
            if x.startswith("("):
                vals.append(float(x.split()[0].lstrip("(")))
            else:
                vals.append(float(x))
        return np.array(vals)
    raise ValueError(f"Could not parse field from {path}")


def main():
    print("=" * 75)
    print("WP 2.4: 1D Planar CJ / ZND Detonation Verification")
    print("Cross-verifying blastFoam with EuDetCode (external/eudetcode)")
    print("=" * 75)

    ref = get_eudetcode_reference()
    D_CJ = ref["D_CJ"]
    P_CJ = ref["P_CJ"]
    T_CJ = ref["T_CJ"]
    rho0 = ref["rho0"]

    print(f"EuDetCode Reference State (COMP-B @ {rho0} kg/m^3):")
    print(f"  Detonation Velocity D_CJ : {D_CJ:.2f} m/s")
    print(f"  Chapman-Jouguet Pressure : {P_CJ/1e9:.3f} GPa ({P_CJ:.3e} Pa)")
    print(f"  Chapman-Jouguet Temp     : {T_CJ:.1f} K")
    print("-" * 75)

    # Run case if time directories not present
    if not os.path.isfile(os.path.join(CASE_DIR, "log.blastFoam")):
        print("Running blastFoam planar_znd simulation...")
        res = subprocess.run(["bash", "-c", f"cd {CASE_DIR} && ./Allrun"], capture_output=True, text=True)
        if res.returncode != 0:
            print(f"Execution failed: {res.stderr[:200]}")
            sys.exit(1)

    log = os.path.join(CASE_DIR, "log.blastFoam")
    if not os.path.isfile(log) or "\nEnd" not in open(log).read():
        print("FAIL: blastFoam did not finish")
        sys.exit(1)

    xc = (np.arange(400) + 0.5) * DX

    # Find written times inside charge domain (< 18 us)
    time_dirs = []
    for d in os.listdir(CASE_DIR):
        p = os.path.join(CASE_DIR, d)
        if os.path.isdir(p):
            try:
                t = float(d)
                if 1e-6 < t <= 1.6e-5:  # Inside solid HE before breakout
                    time_dirs.append((t, d))
            except ValueError:
                pass

    time_dirs.sort(key=lambda x: x[0])

    print(f"{'Time (us)':10} | {'x_num (mm)':12} | {'x_exact (mm)':12} | {'Diff (mm)':10} | {'Rel Diff (%)':12} | Status")
    print("-" * 75)

    all_passed = True
    x_fronts_all = []
    p_front = []

    for t_sec, td in time_dirs:
        lam_path = os.path.join(CASE_DIR, td, "lambda.HE")
        lam = read_openfoam_field(lam_path)

        idx_burned = np.where(lam > 0.5)[0]
        x_num = xc[idx_burned[-1]] if len(idx_burned) > 0 else 0.0
        if len(idx_burned) and 2e-6 <= t_sec <= 1.5e-5:
            p = read_openfoam_field(os.path.join(CASE_DIR, td, "p"))
            i = idx_burned[-1]
            p_front.append(np.max(p[max(0, i - 2):i + 1]))
        x_fronts_all.append(x_num)
        x_exact = D_CJ * t_sec
        diff_mm = (x_num - x_exact) * 1000.0
        rel_diff = abs(x_num - x_exact) / x_exact * 100.0

        # Programmed-burn position: last cell centre reached by the burn
        reached = xc[xc <= x_exact*(1 + 1e-9)]
        x_q = reached[-1] if len(reached) else 0.0
        passed = abs(x_num - x_q) < 1e-9
        if not passed:
            all_passed = False

        status_str = "PASSED" if passed else "FAILED"
        print(f"{t_sec*1e6:10.2f} | {x_num*1000:12.3f} | {x_exact*1000:12.3f} | {diff_mm:10.3f} | {rel_diff:12.3f}% | {status_str}")

    times_arr = [0.0] + [t for t, _ in time_dirs]
    x_arr = [0.0] + [x for x in x_fronts_all]

    slope, intercept = np.polyfit(times_arr, x_arr, 1)
    speed_err = abs(slope - D_CJ) / D_CJ * 100.0
    print(f"Fitted Detonation Velocity     : {slope:.2f} m/s (Ref: {D_CJ:.2f} m/s, Error: {speed_err:.3f}%)")

    D_eos, p_eos, vdet = products_cj_speed()
    eos_err = abs(D_eos/vdet - 1)*100
    print(f"CJ point of the products EOS   : {D_eos:.2f} m/s, {p_eos/1e9:.3f} GPa "
          f"(programmed burn {vdet:.2f} m/s, error {eos_err:.3f}%)")
    eos_ok = eos_err < 0.5
    pf = np.array(p_front)
    front_ok = len(pf) > 0 and np.all(np.abs(pf/P_CJ - 1) < 0.15) and pf.max() <= 1.2*pf.min()
    print(f"Pressure at the burn front     : {pf.min()/1e9:.2f}-{pf.max()/1e9:.2f} GPa "
          f"(p_CJ {P_CJ/1e9:.2f} GPa){'' if front_ok else '  FAILED'}")
    all_passed = all_passed and eos_ok and front_ok

    # Exit check: burn positions, products EOS and front pressure
    if all_passed:
        print("[SUCCESS] WP 2.4 Planar ZND Detonation Verified against EuDetCode.")
        sys.exit(0)
    else:
        print(f"[FAILURE] Planar detonation: speed error {speed_err:.3f}%, products EOS CJ "
              f"error {eos_err:.3f}%, front pressure {'ok' if front_ok else 'out of range'}.")
        sys.exit(1)


if __name__ == "__main__":
    main()

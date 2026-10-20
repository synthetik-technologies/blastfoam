#!/usr/bin/env python3
"""
verify_2d_axi_spherical.py: the spherical Comp-B/Al charge of
tests/1D_spherical/ChargeCompBAl_1D_spherical on a 2D axisymmetric wedge with
AMR (plan step G1, WP 5.1).

Probes on rays at 0, 30, 45, 60 and 90 degrees from the symmetry plane at the
radii of the 1D probes, and every 5 degrees at r = 0.10 and 0.15 m. The rays
on the symmetry plane (0) and on the axis (90 degrees) are reported but not
checked: the cells along each deviate (carbuncle of the HLLC flux at the
grid-aligned shock, felt up to about 10 degrees from them), see UPDATES.md.
Criteria:
1. Angular variation (max - min)/mean in the air of the interior rays (30-60
   degrees at r = 0.08 m; 10-80 degrees at r = 0.10, 0.15 m): peak pressure
   < PEAK_TOL, arrival time < ARRIVAL_TOL.
2. Overlay on the 1D case (its postProcessing must exist, i.e. the 1D test
   has been run): mean of the 30, 45 and 60 degree rays within ONED_TOL_T
   (arrival in the air), ONED_TOL_T_CHARGE (detonation arrival in the charge,
   r = 0.02-0.04 m) and ONED_TOL_P (peak in the air) of the 1D values. The
   tolerances are the change of the 1D case itself from 0.5 to 0.25 mm
   (arrival up to 1.4%, peak up to 2.8%, 2026-10-20).
"""

import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASE_DIR = os.path.join(HERE, "ChargeCompBAl_2D_axi_spherical")
CASE_1D = os.path.join(
    HERE, "..", "1D_spherical", "ChargeCompBAl_1D_spherical"
)
if len(sys.argv) > 1:
    CASE_DIR = sys.argv[1]
if len(sys.argv) > 2:
    CASE_1D = sys.argv[2]

ANGLES = (0, 30, 45, 60, 90)
INTERIOR = (1, 2, 3)            # indices of the 30, 45 and 60 degree rays
RADII = (0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 0.10, 0.15)
CHARGE_RADII = (0.02, 0.03, 0.04)
AIR_RADII = (0.08, 0.10, 0.15)
DENSE_RADII = (0.10, 0.15)
DENSE_ANGLES = tuple(range(0, 91, 5))
P0 = 101325.0

PEAK_TOL = 0.05
ARRIVAL_TOL = 0.02
ONED_TOL_T = 0.015
ONED_TOL_T_CHARGE = 0.015
ONED_TOL_P = 0.05


def read_probes(case, name="pressureProbes"):
    """Time and pressure of all probes (restarts appended in time order)."""
    root = os.path.join(case, "postProcessing", name)
    data = []
    for t in sorted(os.listdir(root), key=float):
        a = np.loadtxt(os.path.join(root, t, "p"), comments="#")
        if data:
            a = a[a[:, 0] > data[-1][-1, 0]]
        data.append(a)
    a = np.vstack(data)
    return a[:, 0], a[:, 1:]


def arrival(t, p, frac=0.5):
    """Time at which p first exceeds p0 + frac*(peak - p0), interpolated."""
    level = P0 + frac*(p.max() - P0)
    i = int(np.argmax(p > level))
    if i == 0:
        return t[0]
    return t[i - 1] + (level - p[i - 1])/(p[i] - p[i - 1])*(t[i] - t[i - 1])


def spread(x):
    return np.ptp(x)/np.mean(x)


def check(ok, good, text):
    print(f"    [{'PASS' if good else 'FAIL'}] {text}")
    return ok and good


def main():
    if not os.path.isdir(os.path.join(CASE_DIR, "postProcessing")):
        print(f"FAIL: no results in {CASE_DIR}")
        return 1

    t2, p2 = read_probes(CASE_DIR)
    nR = len(RADII)
    if p2.shape[1] != nR*len(ANGLES):
        print(f"FAIL: expected {nR*len(ANGLES)} probes, found {p2.shape[1]}")
        return 1

    ta = np.zeros((len(ANGLES), nR))
    pk = np.zeros((len(ANGLES), nR))
    for a in range(len(ANGLES)):
        for j in range(nR):
            p = p2[:, a*nR + j]
            ta[a, j], pk[a, j] = arrival(t2, p), p.max()

    ok = True
    print(f"2D axisymmetric spherical Comp-B/Al (end {t2[-1]*1e6:.1f} us)")
    print("  r [m]  " + "".join(f"{a:>8d}deg" for a in ANGLES))
    for j, r in enumerate(RADII):
        print(f"  {r:5.2f} t " + "".join(f"{v*1e6:11.3f}" for v in ta[:, j])
              + "  [us]")
        print(f"        p " + "".join(f"{v/1e6:11.2f}" for v in pk[:, j])
              + "  [MPa]")

    # 1. Angular variation of the interior rays
    print("  Angular variation (interior rays):")
    j = RADII.index(0.08)
    ok = check(ok, spread(pk[INTERIOR, j]) < PEAK_TOL,
               f"r = 0.08 m, 30-60 deg: peak {spread(pk[INTERIOR, j])*100:.2f}%"
               f" (< {PEAK_TOL*100:g}%)")
    ok = check(ok, spread(ta[INTERIOR, j]) < ARRIVAL_TOL,
               f"r = 0.08 m, 30-60 deg: arrival "
               f"{spread(ta[INTERIOR, j])*100:.2f}% (< {ARRIVAL_TOL*100:g}%)")

    td, pd = read_probes(CASE_DIR, "angularProbes")
    nA = len(DENSE_ANGLES)
    for k, r in enumerate(DENSE_RADII):
        pr = pd[:, k*nA:(k + 1)*nA]
        if pr[:, -2].max() < 2*P0:
            print(f"    r = {r} m: shock not arrived, skipped")
            continue
        pkd = pr.max(axis=0)
        tad = np.array([arrival(td, pr[:, i]) for i in range(nA)])
        inner = slice(2, nA - 2)
        print(f"    r = {r} m: plane {pkd[0]/1e6:.2f} MPa / "
              f"{tad[0]*1e6:.3f} us, axis {pkd[-1]/1e6:.2f} MPa / "
              f"{tad[-1]*1e6:.3f} us, 10-80 deg {pkd[inner].min()/1e6:.2f}-"
              f"{pkd[inner].max()/1e6:.2f} MPa / {tad[inner].min()*1e6:.3f}-"
              f"{tad[inner].max()*1e6:.3f} us (5-85 deg: peak "
              f"{spread(pkd[1:-1])*100:.2f}%; checked: 10-80 deg)")
        ok = check(ok, spread(pkd[inner]) < PEAK_TOL,
                   f"r = {r} m, 10-80 deg: peak {spread(pkd[inner])*100:.2f}%"
                   f" (< {PEAK_TOL*100:g}%)")
        ok = check(ok, spread(tad[inner]) < ARRIVAL_TOL,
                   f"r = {r} m, 10-80 deg: arrival "
                   f"{spread(tad[inner])*100:.2f}% (< {ARRIVAL_TOL*100:g}%)")

    # 2. Overlay on 1D
    if not os.path.isdir(os.path.join(CASE_1D, "postProcessing")):
        print(f"FAIL: no 1D results in {CASE_1D} (run the 1D test first)")
        return 1
    t1, p1 = read_probes(CASE_1D)
    print("  1D overlay (mean of the 30, 45, 60 deg rays vs 1D):")
    for j, r in enumerate(RADII):
        if r not in CHARGE_RADII + AIR_RADII:
            continue
        if r in AIR_RADII and p2[:, INTERIOR[0]*nR + j].max() < 2*P0:
            print(f"    r = {r:4.2f} m: shock not arrived, skipped")
            continue
        t1d, p1d = arrival(t1, p1[:, j]), p1[:, j].max()
        dt = ta[INTERIOR, j].mean()/t1d - 1
        dp = pk[INTERIOR, j].mean()/p1d - 1
        tolT = ONED_TOL_T if r in AIR_RADII else ONED_TOL_T_CHARGE
        good = abs(dt) < tolT and (r not in AIR_RADII or abs(dp) < ONED_TOL_P)
        ok = check(ok, good,
                   f"r = {r:4.2f} m: arrival {ta[INTERIOR, j].mean()*1e6:7.3f}"
                   f" vs {t1d*1e6:7.3f} us ({dt*100:+.2f}%, < {tolT*100:g}%),"
                   f" peak {pk[INTERIOR, j].mean()/1e6:8.2f} vs "
                   f"{p1d/1e6:8.2f} MPa ({dp*100:+.2f}%"
                   + (f", < {ONED_TOL_P*100:g}%)" if r in AIR_RADII else ")"))

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

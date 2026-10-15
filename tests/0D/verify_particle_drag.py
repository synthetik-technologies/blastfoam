#!/usr/bin/env python3
"""
verify_particle_drag.py: Loth (2008) compressible/rarefied drag in
blastEulerFoam (tests/0D/particleDrag_0D, Phase D2).

Particles (Al, alpha 1e-6) start with a velocity U0 in gas at rest; the
slip velocity relaxes through the drag force

    m_p dU_p/dt = -K (U_p - U_g) = -m_g dU_g/dt,
    K = alpha_p 3/4 Cd Re rho_g nu_g/d^2,

with Cd(Re, Ma, Kn) of Loth (2008), implemented here independently and
integrated with RK4. Two runs:
  * compression regime: d = 100 um, 1 bar, U0 = 600 m/s (Ma 1.7 -> 0,
    Re 3900 -> 0);
  * rarefaction regime: d = 1 um, 0.2 bar, U0 = 400 m/s (Re 5, Kn 0.3).
Checks:
1. The correlation has the incompressible limits (Schiller-Naumann with
   the Clift-Gauvin term at Ma -> 0) and continuity checks at Ma = 1.45.
2. Slip velocity history within 1% of U0 of the reference in both runs.
"""

import math
import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.join(HERE, "particleDrag_0D")
WORK = os.path.join(HERE, "particle_drag")

GAMMA = 1.4
CV_GAS = 718.0
MU_GAS = 1.81e-5
RHO0_AL, K0_AL, N_AL, PREF_AL = 2700.0, 72.7e9, 4.3, 101325.0

RUNS = (
    # name, d, p, T, U0, alpha, deltaT, endTime
    ("compression", 1e-4, 1e5, 300.0, 600.0, 1e-6, 1e-7, 2e-3),
    ("rarefaction", 1e-6, 2e4, 300.0, 400.0, 1e-6, 2e-9, 4e-5),
)


def loth_cd(Re, Ma, gamma=GAMMA, TpByT=1.0):
    """Drag coefficient of Loth (2008), AIAA J. 46(9) 2219."""
    if Re < 45.0:
        Kn = math.sqrt(math.pi*gamma/2.0)*Ma/Re
        fKn = 1.0/(1.0 + Kn*(2.514 + 0.8*math.exp(-0.55/Kn))) if Kn > 0 else 1.0
        cd_kn = 24.0/Re*(1.0 + 0.15*Re**0.687)*fKn
        s = max(Ma, 0.01)*math.sqrt(gamma/2.0)
        cd_fm = ((1.0 + 2.0*s**2)*math.exp(-s**2)/(s**3*math.sqrt(math.pi))
                 + (4.0*s**4 + 4.0*s**2 - 1.0)*math.erf(s)/(2.0*s**4)
                 + 2.0/(3.0*s)*math.sqrt(math.pi*TpByT))
        cd_fm_re = cd_fm/(1.0 + (cd_fm/1.63 - 1.0)*math.sqrt(Re/45.0))
        return (cd_kn + Ma**4*cd_fm_re)/(1.0 + Ma**4)
    if Ma < 1.45:
        CM = 5.0/3.0 + 2.0/3.0*math.tanh(3.0*math.log(Ma + 0.1))
    else:
        CM = 2.044 + 0.2*math.exp(-1.8*math.log(Ma/1.5)**2)
    if Ma < 0.89:
        GM = 1.0 - 1.525*Ma**4
    else:
        GM = 0.0002 + 0.0008*math.tanh(12.77*(Ma - 2.02))
    HM = 1.0 - 0.258*CM/(1.0 + 514.0*GM)
    return 24.0/Re*(1.0 + 0.15*Re**0.687)*HM + 0.42*CM/(1.0 + 42500.0*GM/Re**1.16)


def check_limits():
    ok = True
    for Re in (50.0, 300.0, 1e3, 1e4):
        sn = 24.0/Re*(1.0 + 0.15*Re**0.687) + 0.42/(1.0 + 42500.0/Re**1.16)
        cd = loth_cd(Re, 1e-4)
        print(f"  Re {Re:7.0f}, Ma 0: Cd {cd:.4f}, Schiller-Naumann + Clift-Gauvin {sn:.4f}")
        if abs(cd/sn - 1) > 0.01:
            ok = False
    for Re in (100.0, 1e4):
        lo, hi = loth_cd(Re, 1.45 - 1e-9), loth_cd(Re, 1.45 + 1e-9)
        print(f"  Re {Re:7.0f}, Ma 1.45: Cd {lo:.4f} | {hi:.4f}")
        if abs(lo/hi - 1) > 0.02:
            ok = False
    return ok


def env(cmd):
    return ("source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc >/dev/null 2>&1 && "
            "source /opt/OpenFOAM/blastfoam/etc/bashrc >/dev/null 2>&1 && " + cmd)


def run_case(name, d, p, T, U0, alpha, dt, t_end):
    case = os.path.join(WORK, name)
    if os.path.isdir(case):
        shutil.rmtree(case)
    shutil.copytree(CASE, case)
    r = subprocess.run(["bash", "-c", env(f"./Allrun {d!r} {p!r} {T!r} {U0!r} {alpha!r} {dt!r} {t_end!r}")],
                       capture_output=True, text=True, cwd=case)
    logf = os.path.join(case, "log.blastEulerFoam")
    log = open(logf).read() if os.path.exists(logf) else ""
    if r.returncode != 0 or "\nEnd" not in log or "FOAM FATAL" in log:
        print(log[-2000:], r.stdout[-1000:], r.stderr[-1000:])
        sys.exit(f"FAIL: {name}: blastEulerFoam did not finish")
    path = os.path.join(case, "postProcessing", "velocity", "0", "volFieldValue.dat")
    text = open(path).read().replace("(", " ").replace(")", " ")
    data = np.loadtxt(text.splitlines(), comments="#", ndmin=2)
    # columns: t, U.Al (3), U.gas (3), T.Al, T.gas
    return data[:, 0], data[:, 1], data[:, 4], data[:, 8]


def reference(d, p, T, U0, alpha, t_end, n=40000):
    rho_g = p/((GAMMA - 1.0)*CV_GAS*T)
    rho_p = RHO0_AL*(1.0 + N_AL*(p - PREF_AL)/K0_AL)**(1.0/N_AL)
    nu = MU_GAS/rho_g
    c = math.sqrt(GAMMA*(GAMMA - 1.0)*CV_GAS*T)
    m_p, m_g = alpha*rho_p, (1.0 - alpha)*rho_g

    def rhs(y):
        ur = y[0] - y[1]
        Re = max(abs(ur)*d/nu, 1e-12)
        K = alpha*0.75*loth_cd(Re, abs(ur)/c)*Re*rho_g*nu/d**2
        return np.array([-K*ur/m_p, K*ur/m_g])

    y = np.array([U0, 0.0])
    t = np.linspace(0.0, t_end, n + 1)
    ur = np.empty(n + 1)
    dt = t_end/n
    for i in range(n + 1):
        ur[i] = y[0] - y[1]
        if i == n:
            break
        k1 = rhs(y)
        k2 = rhs(y + 0.5*dt*k1)
        k3 = rhs(y + 0.5*dt*k2)
        k4 = rhs(y + dt*k3)
        y = y + dt/6.0*(k1 + 2*k2 + 2*k3 + k4)
    Re0, Ma0 = U0*d/nu, U0/c
    return t, ur, Re0, Ma0, math.sqrt(math.pi*GAMMA/2.0)*Ma0/Re0


def main():
    ok = True
    print("Correlation limits:")
    if not check_limits():
        print("FAIL: Loth correlation limits")
        ok = False

    os.makedirs(WORK, exist_ok=True)
    for name, d, p, T, U0, alpha, dt, t_end in RUNS:
        t, Up, Ug, Tg = run_case(name, d, p, T, U0, alpha, dt, t_end)
        t_ref, ur_ref, Re0, Ma0, Kn0 = reference(d, p, T, U0, alpha, t_end)
        ur = Up - Ug
        err = np.max(np.abs(ur - np.interp(t, t_ref, ur_ref)))/U0
        half = t[np.argmax(ur <= 0.5*U0)] if (ur <= 0.5*U0).any() else float("nan")
        half_ref = t_ref[np.argmax(ur_ref <= 0.5*U0)]
        print(f"{name}: d {d*1e6:g} um, p {p/1e5:g} bar, Re0 {Re0:.1f}, Ma0 {Ma0:.2f}, Kn0 {Kn0:.3g}")
        print(f"  slip halved at {half*1e6:.3f} us (reference {half_ref*1e6:.3f} us), "
              f"final slip {ur[-1]:.3f} m/s (reference {ur_ref[-1]:.3f}), "
              f"max deviation {100*err:.3f}% of U0, gas T change {Tg[-1] - Tg[0]:+.3f} K")
        if err > 0.01:
            print(f"FAIL: {name}: slip history deviates by more than 1% of U0")
            ok = False

    print("SUCCESS: Loth drag matches the reference" if ok else "FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
verify_al_burning.py: Beckstead burning of Al particles in blastEulerFoam
(tests/0D/alBurning_0D, Phase D3).

A dilute Al cloud (alpha 1e-7) at 2500 K (above Tign = 2000 K) burns in air
at 1 atm and 2500 K. The `Beckstead` diameter rate applies the d^n law of
the Beckstead (2005) burn-time correlation

    t_b [ms] = 0.00735 d0^1.8/(X_eff p^0.1 T^0.2),  X_eff = X_O2 (air),

and the burnt mass goes to the gas as Al vapour (no gas chemistry). Runs for
d0 = 1, 10 and 50 um. Checks:
1. Burn time (from the time at which d^n has fallen to 10% of d0^n, i.e.
   0.9 t_b) within 2% of the correlation (exit criterion D3: 5%).
2. d^n decreases linearly in time (deviation < 1% of d0^n).
3. Particle mass follows (d/d0)^3 and the total mass is conserved.
4. Particles below the ignition temperature (1800 K) do not burn.
The cloud is followed until its volume fraction falls below alphaMin
(1e-10, 0.1% of the initial value), where the diameter is reset to d0.
"""

import math
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.join(HERE, "alBurning_0D")
WORK = os.path.join(HERE, "al_burning")

A, N = 0.00735, 1.8
P0, T0, ALPHA = 101325.0, 2500.0, 1e-7
Y_O2, W_O2, W_N2 = 0.2331, 31.9988, 28.0134
D0S = (1e-6, 1e-5, 5e-5)


def burn_time(d0):
    """Beckstead (2005) burn time [s] in air at P0 and T0."""
    X = (Y_O2/W_O2)/(Y_O2/W_O2 + (1 - Y_O2)/W_N2)
    return 1e-3*A*(d0*1e6)**N/(X*(P0/101325.0)**0.1*T0**0.2)


def env(cmd):
    return ("source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc >/dev/null 2>&1 && "
            "source /opt/OpenFOAM/blastfoam/etc/bashrc >/dev/null 2>&1 && " + cmd)


def load(case, name):
    path = os.path.join(case, "postProcessing", name, "0", "volFieldValue.dat")
    return np.loadtxt(path, comments="#", ndmin=2)


def run_case(d0, T=T0, name=None):
    tb = burn_time(d0)
    name = name or f"d{d0*1e6:g}um"
    case = os.path.join(WORK, name)
    if os.path.isdir(case):
        shutil.rmtree(case)
    shutil.copytree(CASE, case)
    dt, t_end = tb/4000, 1.02*tb
    r = subprocess.run(["bash", "-c", env(f"./Allrun {d0!r} {P0!r} {T!r} {ALPHA!r} {dt!r} {t_end!r}")],
                       capture_output=True, text=True, cwd=case)
    logf = os.path.join(case, "log.blastEulerFoam")
    log = open(logf).read() if os.path.exists(logf) else ""
    if r.returncode != 0 or "\nEnd" not in log or "FOAM FATAL" in log:
        print(log[-2000:], r.stdout[-1000:], r.stderr[-1000:])
        sys.exit(f"FAIL: {name}: blastEulerFoam did not finish")
    return name, tb, load(case, "state"), load(case, "mass")


def main():
    ok = True
    os.makedirs(WORK, exist_ok=True)
    for d0 in D0S:
        name, tb, state, mass = run_case(d0)
        t, d = state[:, 0], state[:, 1]
        # Below alphaMin (0.1% of the cloud left) the diameter is reset to
        # d0: use the history before that
        up = np.nonzero(np.diff(d) > 0)[0]
        if len(up):
            t, d, mass = t[:up[0] + 1], d[:up[0] + 1], mass[:up[0] + 1]
        dn = (d/d0)**N
        i90 = np.argmax(dn <= 0.1)
        if dn[i90] > 0.1:
            print(f"FAIL: {name}: the particles did not burn")
            ok = False
            continue
        # interpolate the crossing of dn = 0.1
        t90 = np.interp(0.1, dn[i90::-1][:2], t[i90::-1][:2])
        tb_sim = t90/0.9
        err_tb = tb_sim/tb - 1
        sel = dn > 0.1
        lin = np.max(np.abs(dn[sel] - (1 - t[sel]/tb_sim)))
        # mass: particle mass fraction vs (d/d0)^3 while d > 0.3 d0
        mAl, mg = mass[:, 1], mass[:, 2]
        sel3 = d > 0.3*d0
        m_err = np.max(np.abs(mAl[sel3]/mAl[0] - (d[sel3]/d0)**3))
        m_tot = np.max(np.abs((mAl + mg) - (mAl[0] + mg[0])))/(mAl[0] + mg[0])
        print(f"{name}: t_b {tb_sim*1e3:.5g} ms (Beckstead {tb*1e3:.5g} ms, {100*err_tb:+.2f}%), "
              f"d^n linearity {100*lin:.2f}%, mass vs (d/d0)^3 {100*m_err:.2f}%, "
              f"total mass {m_tot:.1e}, Al burnt {100*(1 - mass[-1, 1]/mass[0, 1]):.2f}%")
        if abs(err_tb) > 0.02:
            print(f"FAIL: {name}: burn time differs by more than 2%")
            ok = False
        if lin > 0.01:
            print(f"FAIL: {name}: d^n not linear in time")
            ok = False
        if m_err > 0.01 or m_tot > 1e-8:
            print(f"FAIL: {name}: particle or total mass inconsistent")
            ok = False

    # Below the ignition temperature (1800 K < Tign) nothing burns
    name, tb, state, mass = run_case(1e-5, T=1800.0, name="d10um_cold")
    burnt = 1 - mass[-1, 1]/mass[0, 1]
    print(f"{name}: particles at 1800 K (< Tign 2000 K): d/d0 {state[-1, 1]/1e-5:.6f}, "
          f"Al burnt {100*burnt:.2e}%")
    if burnt > 1e-10 or abs(state[-1, 1]/1e-5 - 1) > 1e-6:
        print(f"FAIL: {name}: particles below Tign burn")
        ok = False

    print("SUCCESS: Beckstead burning matches the correlation" if ok else "FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

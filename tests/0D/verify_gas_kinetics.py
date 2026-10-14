#!/usr/bin/env python3
"""
verify_gas_kinetics.py: finite-rate gas mechanism below Teq (plan step C6).

The mechanism of the Comp-B/Al gas phase (utilities/mechanism/
make_gas_mechanism.py: GRI-3.0 H2/O2/CO and C1 subset, reversible, with the
nasa9 thermodynamics; local equilibrium above Teq = 2500 K) is run in the
single closed cell tests/0D/gasKinetics_0D (blastFoam, perfect gas) for
mixtures of the Comp-B/Al detonation product gases (pressurePhaseChange
target fractions without C(gr) and Al2O3(s)) with air, and compared with the
full GRI-3.0 mechanism (53 species, 325 reactions) in a Cantera
constant-volume reactor (utilities/mechanism/cantera_cv_ignition.cpp, built
against /usr/local Cantera).

Conditions: 50% and 80% air by mass; 1200, 1500 and 2000 K; 1, 10 and 50 bar.

Pass criteria:
- ignition: the times at which T rises by 50, 200 and 400 K (rises below
  Teq) within 25% of GRI-3.0;
- burnout: final temperature within 1% and the mass fractions of O2, CO,
  CO2 and H2O within 0.01 of GRI-3.0 at the same end time.
"""

import os
import re
import shutil
import subprocess
import sys

import numpy as np

ROOT = "/opt/OpenFOAM/blastfoam"
HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.join(HERE, "gasKinetics_0D")
WORK = os.path.join(HERE, "gas_kinetics")
DRIVER_SRC = os.path.join(ROOT, "utilities/mechanism/cantera_cv_ignition.cpp")
GRI = "/usr/local/share/cantera/data/gri30.yaml"

PRODUCTS = {"H2": 0.00433265, "N2": 0.28739527, "CO": 0.28672579,
            "H2O": 0.11160467, "CO2": 0.15365548, "NH3": 0.00144386,
            "CH4": 0.02801544}
AIR = {"N2": 0.7669, "O2": 0.2331}
RISES = (50, 200, 400)


def env(cmd):
    return (
        "source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc >/dev/null 2>&1 && "
        "source /opt/OpenFOAM/blastfoam/etc/bashrc >/dev/null 2>&1 && " + cmd
    )


def mixture(fair):
    s = sum(PRODUCTS.values())
    Y = {k: (1 - fair)*v/s for k, v in PRODUCTS.items()}
    for k, v in AIR.items():
        Y[k] = Y.get(k, 0) + fair*v
    return Y


def ystr(Y):
    return ", ".join(f"{k}:{v:.10g}" for k, v in Y.items())


def build_driver(exe):
    if os.path.exists(exe) and os.path.getmtime(exe) > os.path.getmtime(DRIVER_SRC):
        return
    subprocess.run(["g++", "-O2", "-std=c++17", "-I/usr/local/include", DRIVER_SRC, "-o", exe,
                    "-L/usr/local/lib", "-lcantera_shared", "-Wl,-rpath,/usr/local/lib"], check=True)


def cantera(exe, T0, p0, Y, tEnd, nOut=4000):
    out = subprocess.run([exe, GRI, repr(T0), repr(p0), ystr(Y), repr(tEnd), str(nOut)],
                         capture_output=True, text=True, check=True).stdout.splitlines()
    hist, Yf = [], {}
    for l in out:
        t = l.split()
        if t[0] == "Y":
            Yf[t[1]] = float(t[2])
        elif t[0] not in ("ignition", "final"):
            hist.append((float(t[0]), float(t[1])))
    h = np.array(hist)
    return h[:, 0], h[:, 1], Yf


def blastfoam(name, T0, p0, Y, dt, tEnd):
    case = os.path.join(WORK, name)
    if os.path.isdir(case):
        shutil.rmtree(case)
    shutil.copytree(CASE, case)
    r = subprocess.run(["bash", "-c", env(f"./Allrun {T0} {p0} '{ystr(Y)}' {dt!r} {tEnd!r}")],
                       capture_output=True, text=True, cwd=case)
    logf = os.path.join(case, "log.blastFoam")
    log = open(logf).read() if os.path.exists(logf) else ""
    if r.returncode != 0 or "\nEnd" not in log or "FOAM FATAL" in log:
        print(log[-2000:], r.stdout[-1000:], r.stderr[-1000:])
        sys.exit(f"FAIL: {name}: blastFoam did not finish")
    t, T, cur = [0.0], [T0], None
    for l in log.splitlines():
        m = re.match(r"^Time = (\S+)", l)
        if m:
            cur = float(m.group(1))
        m = re.match(r"^max\(T\): ([^,\s]+)", l)
        if m and cur is not None:
            t.append(cur)
            T.append(float(m.group(1)))
    last = sorted((d for d in os.listdir(case) if re.fullmatch(r"[0-9.e+-]+", d) and d != "0"),
                  key=float)[-1]
    Yf = {}
    for k in ("O2", "CO", "CO2", "H2O"):
        txt = open(os.path.join(case, last, k)).read()
        Yf[k] = float(re.search(r"internalField\s+uniform\s+([^;]+);", txt).group(1))
    return np.array(t), np.array(T), Yf, "frozen (not converged" in log


def crossing(t, T, Tc):
    w = np.where(T >= Tc)[0]
    if len(w) == 0:
        return float("nan")
    i = w[0]
    return t[i - 1] + (Tc - T[i - 1])/(T[i] - T[i - 1])*(t[i] - t[i - 1])


def main():
    os.makedirs(WORK, exist_ok=True)
    exe = os.path.join(WORK, "cantera_cv_ignition")
    build_driver(exe)

    ok = True
    print(f"{'air':>4s} {'T0':>5s} {'p0 bar':>6s} | " + " ".join(f"{'t+'+str(d)+' GRI':>11s} {'ratio':>5s}" for d in RISES)
          + f" | {'Tf GRI':>7s} {'Tf':>7s} {'dYmax':>6s}")
    for fair in (0.5, 0.8):
        Y = mixture(fair)
        for T0 in (1200.0, 1500.0, 2000.0):
            for p0 in (1e5, 1e6, 5e6):
                # Reference: find the time scale, then run to burnout
                tg, Tg, _ = cantera(exe, T0, p0, Y, 0.5, 4000)
                Tfin = Tg[-1]
                t995 = crossing(tg, Tg, T0 + 0.995*(Tfin - T0))
                tEnd = float(2*t995)
                tg, Tg, Yg = cantera(exe, T0, p0, Y, tEnd, 4000)
                t50 = crossing(tg, Tg, T0 + 50)
                dt = float(min(t50/200, tEnd/4000))
                name = f"air{int(fair*100)}_T{int(T0)}_p{int(p0/1e5)}"
                tb, Tb, Yb, frozen = blastfoam(name, T0, p0, Y, dt, tEnd)

                cells = []
                good = not frozen
                for d in RISES:
                    if T0 + d > 2500 - 50 or T0 + d > Tfin:
                        cells.append(f"{'-':>11s} {'-':>5s}")
                        continue
                    a, b = crossing(tg, Tg, T0 + d), crossing(tb, Tb, T0 + d)
                    rr = b/a
                    good = good and abs(rr - 1) < 0.25
                    cells.append(f"{a:11.3e} {rr:5.2f}")
                dTf = abs(Tb[-1] - Tg[-1])/Tg[-1]
                dY = max(abs(Yb[k] - Yg.get(k, 0)) for k in Yb)
                good = good and dTf < 0.01 and dY < 0.01
                ok = ok and good
                print(f"{fair*100:4.0f} {T0:5.0f} {p0/1e5:6.0f} | " + " ".join(cells)
                      + f" | {Tg[-1]:7.1f} {Tb[-1]:7.1f} {dY:6.4f}{'' if good else '  FAIL'}"
                      + ("  (frozen cells)" if frozen else ""))

    print("SUCCESS: gas mechanism agrees with GRI-3.0" if ok else "FAILURE: gas mechanism vs GRI-3.0")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

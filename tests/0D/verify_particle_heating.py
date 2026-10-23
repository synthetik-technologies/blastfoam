#!/usr/bin/env python3
"""
verify_particle_heating.py: rarefied/compressible particle heat transfer in
blastEulerFoam (tests/0D/particleHeating_0D, Phase D2).

Cold Al particles (alpha 1e-6, 300 K, Cv 900) are heated by gas at 2500 K
through `RanzMarshall` with the `Kavanau` Nusselt number,

    Nu_c = 2 + 0.459 Re^0.55 Pr^0.33,  Nu = Nu_c/(1 + 3.42 Ma Nu_c/(Re Pr)),
    q = alpha_p 6 kappa Nu/d^2 (T_gas - T_p)   [W/m^3],

with kappa = Cp mu/Pr of the gas (const transport, Pr 0.72). Two runs,
integrated independently here (RK4):
  * no slip, d = 1 um at 0.2 bar: Ma/Re = nu/(c d) = 0.65, Nu = 0.28
    (2 in the continuum);
  * slip 300 m/s, d = 10 um at 1 bar: drag (Loth, as in
    verify_particle_drag.py) and heat transfer together, Re 23, Ma 0.3.
Check: particle temperature history within 0.5% of the reference.
"""

import math
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from verify_particle_drag import loth_cd  # noqa: E402

CASE = os.path.join(HERE, "particleHeating_0D")
WORK = os.path.join(HERE, "particle_heating")

GAMMA, CV_GAS, MU_GAS, PR_GAS = 1.4, 718.0, 1.81e-5, 0.72
CV_AL = 900.0
RHO0_AL, K0_AL, N_AL, PREF_AL = 2700.0, 72.7e9, 4.3, 101325.0

RUNS = (
    # name, d, p, T, Tp, U0, alpha, deltaT, endTime
    ("rarefied", 1e-6, 2e4, 2500.0, 300.0, 0.0, 1e-6, 2e-9, 6e-5),
    ("slip", 1e-5, 1e5, 2500.0, 300.0, 300.0, 1e-6, 2e-8, 1.5e-3),
)


def nusselt(Re, Pr, MaByRe):
    Nuc = 2.0 + 0.459*Re**0.55*Pr**0.33
    return Nuc/(1.0 + 3.42*MaByRe*Nuc/Pr)


def env(cmd):
    return ("source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc >/dev/null 2>&1 && "
            "source /opt/OpenFOAM/blastfoam/etc/bashrc >/dev/null 2>&1 && " + cmd)


def run_case(name, d, p, T, Tp, U0, alpha, dt, t_end):
    case = os.path.join(WORK, name)
    if os.path.isdir(case):
        shutil.rmtree(case)
    shutil.copytree(CASE, case)
    r = subprocess.run(["bash", "-c", env(f"./Allrun {d!r} {p!r} {T!r} {Tp!r} {U0!r} {alpha!r} {dt!r} {t_end!r}")],
                       capture_output=True, text=True, cwd=case)
    logf = os.path.join(case, "log.blastEulerFoam")
    log = open(logf).read() if os.path.exists(logf) else ""
    if r.returncode != 0 or "\nEnd" not in log or "FOAM FATAL" in log:
        print(log[-2000:], r.stdout[-1000:], r.stderr[-1000:])
        sys.exit(f"FAIL: {name}: blastEulerFoam did not finish")
    path = os.path.join(case, "postProcessing", "state", "0", "volFieldValue.dat")
    text = open(path).read().replace("(", " ").replace(")", " ")
    data = np.loadtxt(text.splitlines(), comments="#", ndmin=2)
    # columns: t, U.Al (3), U.gas (3), T.Al, T.gas
    return data[:, 0], data[:, 7], data[:, 8]


def reference(d, p, T, Tp, U0, alpha, t_end, n=40000):
    rho_g = p/((GAMMA - 1.0)*CV_GAS*T)
    rho_p = RHO0_AL*(1.0 + N_AL*(p - PREF_AL)/K0_AL)**(1.0/N_AL)
    nu = MU_GAS/rho_g
    kappa = GAMMA*CV_GAS*MU_GAS/PR_GAS
    m_p, m_g = alpha*rho_p, (1.0 - alpha)*rho_g

    def rhs(y):
        Up, Ug, Tp_, Tg = y
        c = math.sqrt(GAMMA*(GAMMA - 1.0)*CV_GAS*Tg)
        ur = Up - Ug
        Re = max(abs(ur)*d/nu, 1e-12)
        K = alpha*0.75*loth_cd(Re, abs(ur)/c, TpByT=Tp_/Tg)*Re*rho_g*nu/d**2 if ur else 0.0
        q = alpha*6.0*kappa*nusselt(Re, PR_GAS, nu/(c*d))/d**2*(Tg - Tp_)
        return np.array([-K*ur/m_p, K*ur/m_g, q/(m_p*CV_AL), -q/(m_g*CV_GAS)])

    y = np.array([U0, 0.0, Tp, T])
    t = np.linspace(0.0, t_end, n + 1)
    Tps = np.empty(n + 1)
    dt = t_end/n
    for i in range(n + 1):
        Tps[i] = y[2]
        if i == n:
            break
        k1 = rhs(y)
        k2 = rhs(y + 0.5*dt*k1)
        k3 = rhs(y + 0.5*dt*k2)
        k4 = rhs(y + dt*k3)
        y = y + dt/6.0*(k1 + 2*k2 + 2*k3 + k4)
    c0 = math.sqrt(GAMMA*(GAMMA - 1.0)*CV_GAS*T)
    return t, Tps, nu/(c0*d), nusselt(U0*d/nu, PR_GAS, nu/(c0*d))


def main():
    ok = True
    os.makedirs(WORK, exist_ok=True)
    for name, d, p, T, Tp, U0, alpha, dt, t_end in RUNS:
        t, T_al, T_g = run_case(name, d, p, T, Tp, U0, alpha, dt, t_end)
        t_ref, T_ref, MaByRe, Nu0 = reference(d, p, T, Tp, U0, alpha, t_end)
        err = np.max(np.abs(T_al - np.interp(t, t_ref, T_ref))/np.interp(t, t_ref, T_ref))
        i_half = np.argmax(T_al >= 0.5*(Tp + T))
        i_half_ref = np.argmax(T_ref >= 0.5*(Tp + T))
        print(f"{name}: d {d*1e6:g} um, p {p/1e5:g} bar, U0 {U0:g} m/s, Ma/Re {MaByRe:.3g}, "
              f"Nu(t=0) {Nu0:.3f}: half-way heating at {t[i_half]*1e6:.3f} us "
              f"(reference {t_ref[i_half_ref]*1e6:.3f} us), max T_Al deviation {100*err:.3f}%")
        if err > 5e-3:
            print(f"FAIL: {name}: particle temperature deviates by more than 0.5%")
            ok = False

    print("SUCCESS: Kavanau heat transfer matches the reference" if ok else "FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

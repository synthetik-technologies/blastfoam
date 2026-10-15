#!/usr/bin/env python3
"""
verify_al_melting.py: Al particle heating through melting in blastEulerFoam
(tests/0D/alMelting_0D, Phase D1).

A closed cell holds hot gas (ideal gas, 50 bar, 2500 K) and 30 um Al
particles (alpha 1e-4, 300 K) with the `eMelting` thermodynamics: solid Cv
below and liquid Cv above Tm = 933.47 K, and the latent heat Lm released over
dTm = 2 K. Heat transfer has a constant coefficient Kh, so the particle
energy per unit volume obeys

    d(m_p e_p)/dt = alpha_Al Kh (T_gas - T_Al) = -d(m_g e_g)/dt,

which is integrated here independently (RK4 on the energies, exact
temperature inversion). Checks:
1. Melting plateau: T_Al stays within the melting range for the time the
   reference takes to supply Lm; start and end of the plateau within 1%.
2. T_Al and T_gas histories within 0.5% of the reference.
3. Total energy of both phases conserved to 1e-8.
"""

import os
import re
import sys

import numpy as np

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
CASE_DIR = os.path.join(TESTS_DIR, "alMelting_0D")
POST_DIR = os.path.join(CASE_DIR, "postProcessing")


def entry(path, key):
    with open(path) as f:
        text = re.sub(r"//.*", "", f.read())
    return float(re.search(rf"\b{key}\s+([-+0-9.eE]+)\s*;", text).group(1))


def internal(name):
    return entry(os.path.join(CASE_DIR, "0", name), "uniform")


def load_columns(name):
    path = os.path.join(POST_DIR, name, "0", "volFieldValue.dat")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing postProcessing file: {path}")
    data = np.loadtxt(path, comments="#", ndmin=2)
    return data[:, 0], data[:, 1:]


class Melting:
    """Thermal energy of eMelting relative to Tref (Esref = Cv Tref)."""

    def __init__(self, path):
        self.Cv = entry(path, "Cv")
        self.CvL = entry(path, "CvLiquid")
        self.Tm = entry(path, "Tm")
        self.Lm = entry(path, "Lm")
        self.dTm = entry(path, "dTm")
        self.T1 = self.Tm - 0.5*self.dTm
        self.T2 = self.Tm + 0.5*self.dTm

    def e(self, T):
        return (self.Cv*min(T, self.T1)
                + (self.Cv + self.Lm/self.dTm)*min(max(T - self.T1, 0.0), self.dTm)
                + self.CvL*max(T - self.T2, 0.0))

    def T(self, e):
        e1 = self.Cv*self.T1
        e2 = e1 + self.Cv*self.dTm + self.Lm
        if e <= e1:
            return e/self.Cv
        if e <= e2:
            return self.T1 + (e - e1)/(self.Cv + self.Lm/self.dTm)
        return self.T2 + (e - e2)/self.CvL


def reference(t_end, n=200000):
    al = Melting(os.path.join(CASE_DIR, "constant", "physicalProperties.Al"))
    cv_g = entry(os.path.join(CASE_DIR, "constant", "physicalProperties.gas"), "Cv")
    with open(os.path.join(CASE_DIR, "constant", "phaseProperties")) as f:
        Kh = float(re.search(r"\bKh\s+([-+0-9.eE]+)\s*;", f.read()).group(1))
    alpha_al = internal("alpha.Al")
    m_p = alpha_al*internal("rho.Al")
    m_g = internal("alpha.gas")*internal("rho.gas")

    def rhs(y):
        q = alpha_al*Kh*(y[1]/(m_g*cv_g) - al.T(y[0]/m_p))
        return np.array([q, -q])

    y = np.array([m_p*al.e(internal("T.Al")), m_g*cv_g*internal("T.gas")])
    dt = t_end/n
    t = np.linspace(0.0, t_end, n + 1)
    T_al = np.empty(n + 1)
    T_g = np.empty(n + 1)
    for i in range(n + 1):
        T_al[i] = al.T(y[0]/m_p)
        T_g[i] = y[1]/(m_g*cv_g)
        if i == n:
            break
        k1 = rhs(y)
        k2 = rhs(y + 0.5*dt*k1)
        k3 = rhs(y + 0.5*dt*k2)
        k4 = rhs(y + dt*k3)
        y = y + dt/6.0*(k1 + 2*k2 + 2*k3 + k4)
    return al, t, T_al, T_g


def plateau(t, T, al):
    """First times at which T reaches the start and the end of melting."""
    i1 = np.argmax(T >= al.T1 + 1e-6)
    i2 = np.argmax(T >= al.T2 - 1e-6)
    return t[i1], t[i2]


def main():
    t, T = load_columns("temperature")
    T_al, T_g = T[:, 0], T[:, 1]
    _, E = load_columns("totalEnergy")
    E = E.sum(axis=1)

    al, t_ref, T_al_ref, T_g_ref = reference(t[-1])
    ok = True

    s, e = plateau(t, T_al, al)
    s_ref, e_ref = plateau(t_ref, T_al_ref, al)
    print(f"Melting start: {s*1e3:.4f} ms (reference {s_ref*1e3:.4f} ms)")
    print(f"Melting end:   {e*1e3:.4f} ms (reference {e_ref*1e3:.4f} ms)")
    for name, v, r in (("start", s, s_ref), ("end", e, e_ref)):
        err = abs(v/r - 1)
        if err > 0.01:
            print(f"FAIL: melting {name} differs by {100*err:.2f}% (> 1%)")
            ok = False
    inside = (t > s) & (t < e)
    if inside.any():
        lo, hi = T_al[inside].min(), T_al[inside].max()
        print(f"T_Al on the plateau: {lo:.2f} - {hi:.2f} K (melting range {al.T1:.2f} - {al.T2:.2f} K)")
        if lo < al.T1 - 0.01 or hi > al.T2 + 0.01:
            print("FAIL: T_Al leaves the melting range during melting")
            ok = False
    else:
        print("FAIL: no melting plateau")
        ok = False

    err_al = np.max(np.abs(T_al - np.interp(t, t_ref, T_al_ref))/np.interp(t, t_ref, T_al_ref))
    err_g = np.max(np.abs(T_g - np.interp(t, t_ref, T_g_ref))/np.interp(t, t_ref, T_g_ref))
    print(f"Max deviation from the reference: T_Al {100*err_al:.3f}%, T_gas {100*err_g:.4f}%")
    print(f"Final T_Al {T_al[-1]:.2f} K (reference {T_al_ref[-1]:.2f}), "
          f"T_gas {T_g[-1]:.2f} K (reference {T_g_ref[-1]:.2f})")
    if err_al > 5e-3 or err_g > 5e-3:
        print("FAIL: temperature history deviates from the reference by more than 0.5%")
        ok = False

    drift = np.max(np.abs(E - E[0]))/abs(E[0])
    print(f"Total energy drift: {drift:.2e}")
    if drift > 1e-8:
        print("FAIL: total energy not conserved (> 1e-8)")
        ok = False

    print("SUCCESS: Al melting plateau matches the reference" if ok else "FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

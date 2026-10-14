#!/usr/bin/env python3
"""
verify_shock_dissociating_air.py: strong shock in dissociating air (Work
Package 2.1, Phase C5).

A high-pressure, high-temperature driver (1700 bar, 7500 K) drives a shock of
about Mach 15 into air at 1 atm and 300 K (tests/1D/shock_dissociating_air,
blastFoam). The air is N2, O2, NO, N and O with NASA-9 thermodynamics to
20 000 K (perfect gas), in local chemical equilibrium above 2500 K
(`equilibrium { Teq 2500; }`, no finite-rate reactions).

Reference: the exact Rankine-Hugoniot jump at the shock speed measured in the
simulation, for the same species data,
- in chemical equilibrium behind the shock (pass criterion: post-shock p, T,
  rho and u within 2.5%), and
- frozen (no reactions; shown for comparison: about 10 000 K).
The equilibrium states of the reference are computed here independently of
the solver (element potentials, Newton on the dual problem). The
post-shock composition of the simulation is compared with the reference.
"""

import os
import re
import subprocess
import sys

import numpy as np
from scipy.optimize import brentq, fsolve

R_U = 8314.462618       # J/(kmol K), OpenFOAM RR
P_STD = 1e5             # standard pressure of the NASA-9 data
CASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "shock_dissociating_air")
SPECIES = ["N2", "O2", "NO", "N", "O"]


# --- Species data (from the case's NASA-9 species file) --------------------

class Species:
    def __init__(self, name, W, Ts, coeffs):
        self.name = name
        self.W = W
        self.Ts = Ts
        self.a = [np.array(c) for c in coeffs]
        self.elements = {}
        for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", name):
            self.elements[e] = self.elements.get(e, 0) + (int(n) if n else 1)

    def _fit(self, T):
        Tf = min(max(T, self.Ts[0]), self.Ts[-1])
        a = self.a[-1]
        for i in range(len(self.a) - 1):
            if Tf < self.Ts[i + 1]:
                a = self.a[i]
                break
        cp = a[0]/Tf**2 + a[1]/Tf + a[2] + Tf*(a[3] + Tf*(a[4] + Tf*(a[5] + Tf*a[6])))
        h = (-a[0]/Tf + a[1]*np.log(Tf) + a[7]
             + Tf*(a[2] + Tf*(a[3]/2 + Tf*(a[4]/3 + Tf*(a[5]/4 + Tf*a[6]/5)))))
        s = (-a[0]/(2*Tf**2) - a[1]/Tf + a[2]*np.log(Tf) + a[8]
             + Tf*(a[3] + Tf*(a[4]/2 + Tf*(a[5]/3 + Tf*a[6]/4))))
        return Tf, cp, h, s

    def cp(self, T):            # J/(kg K)
        return R_U/self.W*self._fit(T)[1]

    def h(self, T):             # absolute enthalpy, J/kg
        Tf, cp, h, _ = self._fit(T)
        return R_U/self.W*(h + cp*(T - Tf))

    def g(self, T):             # G/(R_u T) per kmol
        Tf, cp, h, s = self._fit(T)
        return (h + cp*(T - Tf))/T - (s + cp*np.log(T/Tf))


def air_species(case_dir=CASE_DIR):
    text = open(os.path.join(case_dir, "constant", "thermo.nasa9Air")).read()
    sp = {}
    for m in re.finditer(r'^"?([A-Za-z0-9()]+)"?\s*\n\{(.*?)\n\}', text, re.S | re.M):
        name, body = m.group(1), m.group(2)
        if name not in SPECIES:
            continue
        W = float(re.search(r"molWeight\s+([0-9.eE+-]+);", body).group(1))
        Ts = [float(x) for x in re.search(r"Tswitch\s+\(([^)]*)\)", body).group(1).split()]
        rows = re.search(r"coeffs\s*\((.*?)\);", body, re.S).group(1)
        coeffs = [[float(x) for x in r.split()] for r in re.findall(r"\(([^()]*)\)", rows)]
        sp[name] = Species(name, W, Ts, coeffs)
    return [sp[k] for k in SPECIES]


def air_mass_fractions(case_dir=CASE_DIR):
    Y = {}
    for k in SPECIES:
        m = re.search(r"internalField\s+uniform\s+([^;]+);",
                      open(os.path.join(case_dir, "0_orig", k)).read())
        Y[k] = float(m.group(1))
    return np.array([Y[k] for k in SPECIES])


# --- Equilibrium at fixed T and density (independent of the solver) ---------

def equilibrium_TV(S, b, T, rho, lam=None):
    """Gas amounts [kmol/kg] in equilibrium at T, rho for the element
    abundances b (dict, kmol/kg): Newton on the concave dual."""
    els = sorted(b)
    A = np.array([[s.elements.get(e, 0) for e in els] for s in S], dtype=float)
    bv = np.array([b[e] for e in els])
    g = np.array([s.g(T) for s in S])
    lnc0 = np.log(P_STD/(rho*R_U*T))
    if lam is None:
        lam = np.log(bv) - lnc0
    for _ in range(200):
        n = np.exp(np.minimum(lnc0 - g + A @ lam, 300))
        r = bv - A.T @ n
        H = (A.T*n) @ A
        d = np.linalg.solve(H, r)
        step = min(1.0, 2.0/max(np.max(np.abs(d)), 1e-300))
        lam = lam + step*d
        if step == 1.0 and np.max(np.abs(d)) < 1e-12:
            break
    n = np.exp(lnc0 - g + A @ lam)
    return n, lam


def mixture(S, n, T, rho):
    """p, absolute h [J/kg] and mass fractions for amounts n [kmol/kg]."""
    p = rho*np.sum(n)*R_U*T
    h = sum(n[i]*S[i].W*S[i].h(T) for i in range(len(S)))
    return p, h, np.array([n[i]*S[i].W for i in range(len(S))])


def abundances(S, Y):
    b = {}
    for s, y in zip(S, Y):
        for e, c in s.elements.items():
            b[e] = b.get(e, 0) + c*y/s.W
    return b


def frozen_state(S, Y, T, rho):
    n = np.array([y/s.W for s, y in zip(S, Y)])
    return mixture(S, n, T, rho)


def jump(Us, p1, T1, S, Y1, equilibrium=True):
    """Exact Rankine-Hugoniot jump at shock speed Us into the gas (p1, T1, Y1)
    at rest: frozen composition or chemical equilibrium behind the shock."""
    n1 = np.array([y/s.W for s, y in zip(S, Y1)])
    rho1 = p1/(np.sum(n1)*R_U*T1)
    h1 = mixture(S, n1, T1, rho1)[1]
    b = abundances(S, Y1)
    lam = [None]

    def state(T, rho):
        if equilibrium:
            n, lam[0] = equilibrium_TV(S, b, T, rho, lam[0])
            return mixture(S, n, T, rho)
        return frozen_state(S, Y1, T, rho)

    def residual(x):
        T, rho = np.exp(x)
        p, h, _ = state(T, rho)
        pRH = p1 + rho1*Us**2*(1 - rho1/rho)
        hRH = h1 + 0.5*Us**2*(1 - (rho1/rho)**2)
        return [p/pRH - 1, (h - hRH)/(0.5*Us**2)]

    # Start from the frozen perfect-gas (gamma 1.4) jump
    a1 = np.sqrt(1.4*np.sum(n1)*R_U*T1)
    M = Us/a1
    rr = 2.4*M**2/(0.4*M**2 + 2)
    T0 = T1*(1 + 2.8/2.4*(M**2 - 1))/rr
    x = fsolve(residual, np.log([T0, rr*rho1]), xtol=1e-13)
    T2, rho2 = np.exp(x)
    p2, h2, Y2 = state(T2, rho2)
    err = np.max(np.abs(residual(x)))
    return {"Us": Us, "rho1": rho1, "p2": p2, "T2": T2, "rho2": rho2,
            "u2": Us*(1 - rho1/rho2), "Y": dict(zip(SPECIES, Y2)), "residual": err}


def sound_speed(S, Y, T, p):
    n = np.array([y/s.W for s, y in zip(S, Y)])
    R = np.sum(n)*R_U
    cp = sum(y*s.cp(T) for s, y in zip(S, Y))
    return np.sqrt(cp/(cp - R)*R*T)


# --- Simulation -------------------------------------------------------------

def read_openfoam_field(file_path, n_cells=400):
    with open(file_path) as f:
        lines = [line.strip() for line in f if line.strip()]
    for line in lines:
        if line.startswith("internalField"):
            parts = line.split()
            if len(parts) >= 3 and parts[1] == "uniform":
                return np.full(n_cells, float(parts[2].rstrip(";").lstrip("(").split()[0]))
    start = lines.index("(") + 1
    end = lines.index(")")
    return np.array([float(x.split()[0].lstrip("(")) for x in lines[start:end]])


def run_and_extract_blastfoam_case():
    def time_dirs():
        out = []
        for entry in os.listdir(CASE_DIR):
            try:
                if os.path.isdir(os.path.join(CASE_DIR, entry)) and float(entry) > 0:
                    out.append((float(entry), entry))
            except ValueError:
                pass
        return sorted(out)

    if not time_dirs():
        print("[INFO] Running blastFoam shock_dissociating_air simulation...")
        res = subprocess.run(["bash", "-c", f"cd {CASE_DIR} && ./Allrun"], capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"blastFoam run failed: {res.stderr}")
    tds = time_dirs()

    def shock_state(t_dir):
        path = os.path.join(CASE_DIR, t_dir)
        fields = {f: read_openfoam_field(os.path.join(path, f)) for f in ("p", "T", "rho", "U") + tuple(SPECIES)}
        p = fields["p"]
        n = len(p)
        x = (np.arange(n) + 0.5)/n
        # Shock: last cell above 1.5 bar; its position is where p crosses the
        # mid-jump value (linear interpolation)
        i_s = np.where(p > 1.5e5)[0][-1]
        p_mid = 0.5*(p[-1] + np.median(p[max(0, i_s - 12):i_s - 5]))
        i = np.where(p > p_mid)[0][-1]
        x_s = x[i] + (p[i] - p_mid)/(p[i] - p[i + 1])*(x[i + 1] - x[i])
        return fields, x, i_s, x_s

    latest_t, latest_dir = tds[-1]
    fields, xc, shock_cell, x_shock = shock_state(latest_dir)

    # Shock speed: least-squares fit of the shock position over the last three writes
    fit = tds[-3:]
    Us = np.polyfit([t for t, _ in fit], [shock_state(d)[3] for _, d in fit], 1)[0]

    # Post-shock state: median of cells 6-12 behind the shock front
    plateau = slice(max(0, shock_cell - 12), shock_cell - 5)
    post = {f: np.median(fields[f][plateau]) for f in fields}
    out = {
        "t": latest_t, "x": xc, "x_shock": x_shock, "Us": Us,
        "p1": fields["p"][-1], "T1": fields["T"][-1], "rho1": fields["rho"][-1],
        "p2": post["p"], "T2": post["T"], "rho2": post["rho"], "u2": post["U"],
        "Y2": {k: post[k] for k in SPECIES}, "fields": fields,
    }
    return out


def main():
    print("=" * 85)
    print("WP 2.1 / Phase C5: strong shock in dissociating air")
    print("blastFoam (NASA-9 air, local equilibrium above 2500 K) vs exact Rankine-Hugoniot jumps")
    print("=" * 85)

    cfd = run_and_extract_blastfoam_case()
    log = os.path.join(CASE_DIR, "log.blastFoam")
    if not os.path.isfile(log) or "\nEnd" not in open(log).read():
        print("FAIL: blastFoam did not finish")
        sys.exit(1)
    frozenCells = "frozen (not converged" in open(log).read()

    S = air_species()
    Y1 = air_mass_fractions()
    Ma = cfd["Us"]/sound_speed(S, Y1, cfd["T1"], cfd["p1"])

    print(f"  Simulation time        : {cfd['t']*1e6:.1f} us")
    print(f"  Shock position         : {cfd['x_shock']:.4f} m")
    print(f"  Shock speed            : {cfd['Us']:.1f} m/s (Ma {Ma:.2f})")
    print(f"  Post-shock p, T, rho, u: {cfd['p2']/1e5:.2f} bar, {cfd['T2']:.1f} K, "
          f"{cfd['rho2']:.4f} kg/m3, {cfd['u2']:.1f} m/s")

    eq = jump(cfd["Us"], cfd["p1"], cfd["T1"], S, Y1, equilibrium=True)
    fr = jump(cfd["Us"], cfd["p1"], cfd["T1"], S, Y1, equilibrium=False)

    print(f"\n  {'Quantity':14} | {'blastFoam':10} | {'equilibrium':11} | {'error':7} | {'frozen':10}")
    print("  " + "-"*64)
    ok = True
    for key, name, scale in (("p2", "p2 [bar]", 1e-5), ("T2", "T2 [K]", 1.0),
                             ("rho2", "rho2 [kg/m3]", 1.0), ("u2", "u2 [m/s]", 1.0)):
        err = abs(cfd[key] - eq[key])/abs(eq[key])*100
        ok = ok and err < 2.5
        print(f"  {name:14} | {cfd[key]*scale:10.4g} | {eq[key]*scale:11.4g} | {err:6.2f}% | {fr[key]*scale:10.4g}")
    print(f"  Composition behind the shock (mass fractions), blastFoam / equilibrium jump:")
    print("    " + ", ".join(f"{k} {cfd['Y2'][k]:.4f}/{eq['Y'][k]:.4f}" for k in SPECIES))
    dY = max(abs(cfd["Y2"][k] - eq["Y"][k]) for k in SPECIES)
    ok = ok and dY < 0.02 and eq["residual"] < 1e-9 and not frozenCells
    print(f"  max |dY| {dY:.4f} (tolerance 0.02); reference residual {eq['residual']:.1e}"
          + ("; FAIL: frozen cells in the log" if frozenCells else ""))

    # Reference over a range of Mach numbers (no CFD)
    print("\n  Equilibrium and frozen jumps into air at 1 atm, 300 K (same species data):")
    print(f"  {'Ma':5} | {'Us (m/s)':9} | {'T2_eq (K)':9} | {'T2_fr (K)':9} | {'rho2/rho1 eq':12} | "
          f"{'Y_O':6} | {'Y_N':6} | {'Y_NO':6}")
    a1 = sound_speed(S, Y1, 300.0, 101325.0)
    for M in (3.0, 5.0, 8.0, 10.0, 15.0):
        e = jump(M*a1, 101325.0, 300.0, S, Y1, True)
        f = jump(M*a1, 101325.0, 300.0, S, Y1, False)
        ok = ok and e["residual"] < 1e-9
        print(f"  {M:5.1f} | {M*a1:9.1f} | {e['T2']:9.1f} | {f['T2']:9.1f} | {e['rho2']/e['rho1']:12.2f} | "
              f"{e['Y']['O']:6.3f} | {e['Y']['N']:6.3f} | {e['Y']['NO']:6.3f}")

    flds = cfd["fields"]
    np.savez(os.path.join(CASE_DIR, "shock_fields.npz"),
             x=cfd["x"], p=flds["p"], T=flds["T"], rho=flds["rho"], U=flds["U"],
             t=cfd["t"], x_shock=cfd["x_shock"], Us=cfd["Us"], Ma_s=Ma,
             p2=cfd["p2"], T2=cfd["T2"], rho2=cfd["rho2"], u2=cfd["u2"],
             p1=cfd["p1"], T1=cfd["T1"], rho1=cfd["rho1"],
             T2_frozen=fr["T2"], T2_eq=eq["T2"], rho2_frozen=fr["rho2"], rho2_eq=eq["rho2"],
             **{f"Y_{k}": flds[k] for k in SPECIES})

    print("=" * 85)
    if ok:
        print("\n[SUCCESS] WP 2.1: post-shock state matches the equilibrium jump conditions.")
        sys.exit(0)
    print("\n[FAILURE] WP 2.1: post-shock state does not match the equilibrium jump conditions.")
    sys.exit(1)


if __name__ == "__main__":
    main()

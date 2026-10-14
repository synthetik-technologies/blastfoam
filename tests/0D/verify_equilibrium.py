#!/usr/bin/env python3
"""
verify_equilibrium.py: chemical equilibrium at fixed density and internal
energy (plan step C3).

The blastFoam solver (blastChemicalEquilibrium, through applications/test/
Test-blastEquilibrium) is compared with

- the EuDetCode equilibrium solver (ideal gas) for gas-only species sets,
  with a gas database written from the same NASA-9 coefficients as the
  blastFoam species file (standard pressure 1e5 Pa); the driver
  utilities/thermo_fitter/eudet_equilibrium/eudet_uv_equilibrium.cpp is
  compiled against external/eudetcode/build. The element abundances are
  passed directly (no ingredient lookup, no default inert Al). States where
  the EuDetCode UV iteration does not converge (its pressure iteration is
  bounded below by 1e5 Pa) have no EuDetCode reference ("-");
- an independent Python solution for all states (concave dual maximised
  with SciPy SLSQP, temperature by Brent's method), which also covers the
  condensed species C(gr) and Al2O3(s) (EuDetCode models them with its own
  condensed EOS).

Cases: air (N2 O2 NO N O); the Comp-B/Al detonation products (gas part of
the pressurePhaseChange target fractions, with 1% Al vapour) alone and mixed
with air (gas species O2 N2 CO H2O CO2 NH3 CH4 C2H6 Al H2 O N NO H OH AlO);
and the full products with C(gr) and Al2O3(s). Densities 0.01-100 kg/m^3,
equilibrium temperatures 2500-20 000 K.

Pass criteria (plan: within 0.5%):
- temperature and pressure within 0.5%, and the energy of the EuDetCode
  solution within 0.5% of the target;
- mass fractions within 0.5% (relative) for Y > 1e-3 and within 1e-5
  (absolute) otherwise;
- blastFoam solution: converged, element balance within 1e-12 and energy
  within 1e-8 (relative to max(|e|, 1 MJ/kg)).
"""

import os
import re
import subprocess
import sys

import numpy as np

ROOT = "/opt/OpenFOAM/blastfoam"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "equilibrium")
EUDET = os.path.join(ROOT, "external/eudetcode")
DRIVER_SRC = os.path.join(ROOT, "utilities/thermo_fitter/eudet_equilibrium/eudet_uv_equilibrium.cpp")
BASE = os.path.join(HERE, "../1D_spherical/ChargeCompBAl_1D_spherical_reacting/constant/thermo.compressibleGas")

RU = 8314.462618      # OpenFOAM RR [J/kmol/K]
PSTD = 1e5

AIR_SPECIES = ["N2", "O2", "NO", "N", "O"]
GAS_SPECIES = ["O2", "N2", "CO", "H2O", "CO2", "NH3", "CH4", "C2H6", "Al",
               "H2", "O", "N", "NO", "H", "OH", "AlO"]
COND_SPECIES = GAS_SPECIES + ["C(gr)", "Al2O3(s)"]
AIR = {"N2": 0.7669, "O2": 0.2331}
PRODUCTS = {"H2": 0.00433265, "N2": 0.28739527, "CO": 0.28672579,
            "H2O": 0.11160467, "CO2": 0.15365548, "NH3": 0.00144386,
            "CH4": 0.02801544, "C(gr)": 0.04906294, "Al2O3(s)": 0.07776391}


def env(cmd):
    return (
        "source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc >/dev/null 2>&1 && "
        "source /opt/OpenFOAM/blastfoam/etc/bashrc >/dev/null 2>&1 && " + cmd
    )


def normalise(Y):
    s = sum(Y.values())
    return {k: v/s for k, v in Y.items()}


def blend(a, b, fa):
    out = {}
    for k, v in a.items():
        out[k] = out.get(k, 0) + fa*v
    for k, v in b.items():
        out[k] = out.get(k, 0) + (1 - fa)*v
    return out


def gas_products():
    Y = {k: v for k, v in PRODUCTS.items() if "(" not in k}
    Y = normalise(Y)
    Y = {k: 0.99*v for k, v in Y.items()}
    Y["Al"] = 0.01
    return Y


# --- Species data --------------------------------------------------------

def read_species(path):
    """NASA-9 coefficients (per kmol, /R), ranges and molar masses."""
    text = open(path).read()
    data = {}
    for m in re.finditer(r'^"?([A-Za-z0-9()]+)"?\s*\n\{(.*?)\n\}', text, re.S | re.M):
        name, body = m.group(1), m.group(2)
        if name == "species":
            continue
        W = float(re.search(r"molWeight\s+([0-9.eE+-]+);", body).group(1))
        Ts = [float(x) for x in re.search(r"Tswitch\s+\(([^)]*)\)", body).group(1).split()]
        rows = re.search(r"coeffs\s*\((.*?)\);", body, re.S).group(1)
        coeffs = [[float(x) for x in r.split()] for r in re.findall(r"\(([^()]*)\)", rows)]
        data[name] = {"W": W, "Tswitch": Ts, "coeffs": coeffs}
    return data


def elements(name):
    comp = {}
    for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", name.split("(")[0]):
        comp[e] = comp.get(e, 0) + (int(n) if n else 1)
    return comp


class Species:
    def __init__(self, name, d):
        self.name = name
        self.W = d["W"]
        self.Ts = d["Tswitch"]
        self.a = [np.array(c) for c in d["coeffs"]]
        self.condensed = "(" in name
        self.elements = elements(name)

    def _rng(self, T):
        for i in range(len(self.a) - 1):
            if T < self.Ts[i + 1]:
                return self.a[i]
        return self.a[-1]

    def _fit(self, T):
        Tf = min(max(T, self.Ts[0]), self.Ts[-1])
        a = self._rng(Tf)
        cp = a[0]/Tf**2 + a[1]/Tf + a[2] + a[3]*Tf + a[4]*Tf**2 + a[5]*Tf**3 + a[6]*Tf**4
        h = (-a[0]/Tf + a[1]*np.log(Tf) + a[7] + a[2]*Tf + a[3]*Tf**2/2
             + a[4]*Tf**3/3 + a[5]*Tf**4/4 + a[6]*Tf**5/5)
        s = (-a[0]/(2*Tf**2) - a[1]/Tf + a[2]*np.log(Tf) + a[8] + a[3]*Tf
             + a[4]*Tf**2/2 + a[5]*Tf**3/3 + a[6]*Tf**4/4)
        return Tf, cp, h, s

    def h(self, T):          # H/R [K]
        Tf, cp, h, _ = self._fit(T)
        return h + cp*(T - Tf)

    def g(self, T):          # G/(R T)
        Tf, cp, h, s = self._fit(T)
        return (h + cp*(T - Tf))/T - (s + cp*np.log(T/Tf))

    def e(self, T):          # absolute internal energy [J/kg], blastFoam model
        return RU/self.W*(self.h(T) - T)


# --- References ----------------------------------------------------------

def write_eudet_db(sp, names, path):
    import json
    db = {"species": {}, "p0": PSTD}
    for n in names:
        s = sp[n]
        db["species"][n] = {"elements": s.elements,
                            "nasaParams": {"coeffs": [list(map(float, a)) for a in s.a],
                                           "Tswitch": s.Ts}}
    json.dump(db, open(path, "w"))
    json.dump({}, open(os.path.join(os.path.dirname(path), "condensed_empty.json"), "w"))


def build_driver(exe):
    if os.path.exists(exe) and os.path.getmtime(exe) > os.path.getmtime(DRIVER_SRC):
        return
    b = os.path.join(EUDET, "build")
    inc = " ".join(f"-I{EUDET}/src/{d}" for d in ["", "ThermoChemistry", "ThermoChemistry/EosMixture",
                                                    "ThermoChemistry/EosPure", "ThermoChemistry/Kinetics",
                                                    "Numerics"])
    libs = " ".join([
        f"-Wl,--whole-archive {b}/src/ThermoChemistry/libeudetcode_thermochemistry.a -Wl,--no-whole-archive",
        "/usr/lib64/libnlopt.so.1.0.0", f"{b}/src/Numerics/libeudetcode_numerics.a",
        "/usr/lib64/libboost_json.so.1.90.0", "/usr/lib64/libboost_program_options.so.1.90.0",
        "/usr/lib64/libboost_container.so.1.90.0"])
    cmd = (f"c++ -O2 -std=gnu++20 -DNDEBUG -DUSE_NLOPT {inc} -isystem /usr/include/eigen3 "
           f"{DRIVER_SRC} -o {exe} {libs}")
    subprocess.run(cmd, shell=True, check=True)


def parse_states(text, nfield, tags=("state",)):
    states = []
    cur = None
    for line in text.splitlines():
        t = line.split()
        if not t:
            continue
        if t[0] in tags:
            cur = {"tag": t[0],
                   "vals": [float(x) if x not in ("true", "false") else float(x == "true")
                            for x in t[1:1 + nfield]], "Y": {}}
            states.append(cur)
        elif cur is not None and len(t) == 2:
            try:
                cur["Y"][t[0]] = float(t[1])
            except ValueError:
                pass
    return states


def eudet(exe, gasdb, sp, Y0, states):
    """EuDetCode UV equilibrium at (rho, e) for the element abundances of
    Y0; returns T, p, Y, converged and the energy of its solution."""
    b = {}
    for k, y in Y0.items():
        for e, c in sp[k].elements.items():
            b[e] = b.get(e, 0) + 1e3*c*y/sp[k].W     # mol/kg
    els = " ".join(f"{e}:{v:.17g}" for e, v in b.items() if v > 0)
    args = [exe, gasdb, os.path.join(os.path.dirname(gasdb), "condensed_empty.json"), els]
    for rho, e in states:
        args += [repr(rho), repr(e)]
    out = subprocess.run(args, capture_output=True, text=True, cwd=EUDET)
    if out.returncode != 0:
        print(out.stdout[-1000:], out.stderr[-2000:])
        sys.exit("FAIL: EuDetCode driver failed")
    return [{"T": r["vals"][2], "p": r["vals"][3], "Y": r["Y"], "conv": r["vals"][4] == 1,
             "e": r["vals"][5]} for r in parse_states(out.stdout, 6)]


def python_uv(sp, names, Y0, rho, ea, T0):
    """Independent equilibrium: concave dual with SLSQP, T by Brent."""
    from scipy.optimize import minimize, brentq
    S = [sp[n] for n in names]
    els = sorted({e for s in S for e in s.elements})
    A = np.array([[s.elements.get(e, 0) for e in els] for s in S], dtype=float)
    n0 = np.array([Y0.get(s.name, 0)/s.W for s in S])
    b = A.T @ n0
    gas = np.array([not s.condensed for s in S])
    V = 1/rho

    def tv(T):
        g = np.array([s.g(T) for s in S])
        lnc0 = np.log(V*PSTD/(RU*T))
        scale = b.sum()

        def f(lam):
            x = np.exp(np.minimum(lnc0 - g[gas] + A[gas] @ lam, 200))
            return (x.sum() - b @ lam)/scale, (A[gas].T @ x - b)/scale

        cons = [{"type": "ineq", "fun": (lambda lam, k=k: g[k] - A[k] @ lam),
                 "jac": (lambda lam, k=k: -A[k])} for k in np.where(~gas)[0]]
        # start from a gas-phase least-squares fit
        nn = np.maximum(n0, 0) + 1e-10*scale
        t = np.log(nn[gas]) - lnc0 + g[gas]
        w = nn[gas]
        lam0 = np.linalg.lstsq(A[gas]*np.sqrt(w)[:, None], t*np.sqrt(w), rcond=None)[0]
        best = None
        for _ in range(3):
            r = minimize(f, lam0, jac=True, constraints=cons, method="SLSQP",
                         options={"ftol": 1e-15, "maxiter": 2000})
            lam0 = r.x
            best = r
        lam = best.x
        n = np.zeros(len(S))
        n[gas] = np.exp(lnc0 - g[gas] + A[gas] @ lam)
        res = b - A[gas].T @ n[gas]
        cond = np.where(~gas)[0]
        if len(cond):
            act = [k for k in cond if g[k] - A[k] @ lam < 1e-6]
            if act:
                nc = np.linalg.lstsq(A[act].T, res, rcond=None)[0]
                n[act] = np.maximum(nc, 0)
        return n

    def energy(T):
        n = tv(T)
        return sum(n[i]*S[i].W*S[i].e(T) for i in range(len(S))) - ea

    T = brentq(energy, max(300, 0.3*T0), 3*T0 + 1000, xtol=1e-9, rtol=1e-12)
    n = tv(T)
    Y = {S[i].name: n[i]*S[i].W for i in range(len(S))}
    p = sum(n[i] for i in range(len(S)) if not S[i].condensed)*RU*T/V
    return T, p, Y


# --- blastFoam -----------------------------------------------------------

def blastfoam(species_file, names, Y0, rhos, Tfs):
    cmd = (f"Test-blastEquilibrium -file {species_file} -species '({' '.join(names)})' "
           f"-mixture '({' '.join(f'({k} {v:.17g})' for k, v in Y0.items() if v > 0)})' "
           f"-rho '({' '.join(repr(r) for r in rhos)})' -Tfrozen '({' '.join(repr(t) for t in Tfs)})'")
    out = subprocess.run(["bash", "-c", env(cmd)], capture_output=True, text=True, cwd=WORK)
    if out.returncode != 0:
        print(out.stdout[-2000:], out.stderr[-2000:])
        sys.exit("FAIL: Test-blastEquilibrium failed")
    return parse_states(out.stdout, 8)


def deviations(s, Tr, pr, Yr, names):
    T, p = s["vals"][2], s["vals"][3]
    dT = abs(T - Tr)/Tr
    dp = abs(p - pr)/pr
    dmaj = max([abs(s["Y"][n] - Yr.get(n, 0))/Yr[n] for n in names if Yr.get(n, 0) > 1e-3] or [0])
    dmin = max([abs(s["Y"][n] - Yr.get(n, 0)) for n in names if Yr.get(n, 0) <= 1e-3] or [0])
    return dT, dp, dmaj, dmin


def compare(label, bf, ref, names):
    """ref: UV solution (T, p, Y, conv) and optionally the TV solution at the
    blastFoam temperature (p, e, Y), used where the reference UV solve did
    not converge."""
    worst = {"T": 0, "p": 0, "Ymaj": 0, "Ymin": 0, "e": 0}
    ok = True
    nuv = 0
    print(f"\n{label}")
    print(f"  {'rho':>8s} {'T':>9s} {'ref':>3s} {'dT/T':>8s} {'dp/p':>8s} {'dY maj':>8s} {'dY min':>8s} "
          f"{'de/e':>8s} {'elem':>7s} {'energy':>7s} {'cpu s':>7s}")
    for s, r in zip(bf, ref):
        rho, ea, T, p, conv, bErr, eErr, cpu = s["vals"]
        scale = max(abs(ea), 1e6)
        if r.get("conv", True):
            kind = "UV"
            dT, dp, dmaj, dmin = deviations(s, r["T"], r["p"], r["Y"], names)
            # energy of the reference solution (its convergence)
            de = abs(r.get("e", ea) - ea)/scale
            nuv += 1
        else:
            kind = "-"
            dT = dp = dmaj = dmin = de = 0
        good = (conv == 1 and dT < 5e-3 and dp < 5e-3 and dmaj < 5e-3 and dmin < 1e-5
                and de < 5e-3 and bErr < 1e-12 and eErr < 1e-8*scale)
        ok = ok and good
        for k, v in zip(("T", "p", "Ymaj", "Ymin", "e"), (dT, dp, dmaj, dmin, de)):
            worst[k] = max(worst[k], v)
        print(f"  {rho:8.3g} {T:9.1f} {kind:>3s} {dT:8.1e} {dp:8.1e} {dmaj:8.1e} {dmin:8.1e} "
              f"{de:8.1e} {bErr:7.0e} {eErr/scale:7.0e} {cpu:7.4f}{'' if good else '  FAIL'}")
    return ok, worst, nuv


def main():
    os.makedirs(WORK, exist_ok=True)
    species_file = os.path.join(WORK, "thermo.nasa9Gas")
    subprocess.run([sys.executable, os.path.join(ROOT, "utilities/thermo_fitter/nasa9_species.py"),
                    "--base", BASE, "--species", *COND_SPECIES, "-o", species_file],
                   check=True, capture_output=True)
    sp = {k: Species(k, v) for k, v in read_species(species_file).items()}

    app = os.path.join(ROOT, "applications/test/Test-blastEquilibrium")
    r = subprocess.run(["bash", "-c", env("which Test-blastEquilibrium")], capture_output=True, text=True)
    if r.returncode != 0:
        subprocess.run(["bash", "-c", env(f"wmake {app}")], check=True, capture_output=True)

    exe = os.path.join(WORK, "eudet_uv_equilibrium")
    build_driver(exe)

    rhos = [0.01, 0.1, 1.0, 10.0, 100.0]
    ok = True
    summary = []

    gas_cases = [
        ("air", AIR_SPECIES, AIR, [3500, 6000, 12000, 25000, 50000, 80000]),
        ("products + 1% Al", GAS_SPECIES, gas_products(), [3000, 5000, 10000, 20000, 40000]),
        ("50% products / 50% air", GAS_SPECIES, blend(gas_products(), AIR, 0.5), [3000, 6000, 12000, 30000]),
        ("20% products / 80% air", GAS_SPECIES, blend(gas_products(), AIR, 0.2), [3500, 8000, 20000, 50000]),
    ]
    for label, names, Y0, Tfs in gas_cases:
        Y0 = normalise(Y0)
        states = [(r_, t) for r_ in rhos for t in Tfs]
        bf = blastfoam(species_file, names, Y0, [s[0] for s in states], [s[1] for s in states])
        # Equilibrium temperatures in the range of the data and above Teq
        keep = [i for i, s in enumerate(bf) if 2500 <= s["vals"][2] <= 20000]
        bf = [bf[i] for i in keep]
        gasdb = os.path.join(WORK, "gas_" + re.sub(r"\W+", "_", label) + ".json")
        write_eudet_db(sp, names, gasdb)
        refs = eudet(exe, gasdb, sp, Y0, [(s["vals"][0], s["vals"][1]) for s in bf])
        good, worst, nuv = compare(f"{label}: vs EuDetCode", bf, refs, names)
        ok = ok and good
        summary.append((label, "EuDetCode", len(bf), nuv, worst, good))
        pyrefs = []
        for s in bf:
            T, p, Y = python_uv(sp, names, Y0, s["vals"][0], s["vals"][1], s["vals"][2])
            pyrefs.append({"T": T, "p": p, "Y": Y})
        good, worst, nuv = compare(f"{label}: vs Python", bf, pyrefs, names)
        ok = ok and good
        summary.append((label, "Python", len(bf), nuv, worst, good))

    # Condensed species: products with C(gr) and Al2O3(s)
    Y0 = normalise(dict(PRODUCTS))
    Tfs = [2500, 3000, 4000, 6000]
    states = [(r_, t) for r_ in [0.1, 1.0, 10.0, 100.0] for t in Tfs]
    bf = blastfoam(species_file, COND_SPECIES, Y0, [s[0] for s in states], [s[1] for s in states])
    bf = [s for s in bf if 2500 <= s["vals"][2] <= 20000]
    refs = []
    for s in bf:
        T, p, Y = python_uv(sp, COND_SPECIES, Y0, s["vals"][0], s["vals"][1], s["vals"][2])
        refs.append({"T": T, "p": p, "Y": Y})
    good, worst, nuv = compare("products with C(gr), Al2O3(s): vs Python (SLSQP)", bf, refs, COND_SPECIES)
    for s in bf:
        cond = {n: s["Y"][n] for n in ("C(gr)", "Al2O3(s)")}
        print(f"  rho {s['vals'][0]:g} T {s['vals'][2]:.0f}: "
              + ", ".join(f"{k} {v:.4f}" for k, v in cond.items()))
    ok = ok and good
    summary.append(("products + condensed", "Python", len(bf), nuv, worst, good))

    print("\nSummary (max deviations)")
    for label, ref, n, nuv, w, good in summary:
        print(f"  {label:24s} vs {ref:9s} {n:3d} states ({nuv:2d} UV): T {w['T']:.1e}, p {w['p']:.1e}, "
              f"Y(>1e-3) {w['Ymaj']:.1e}, Y(<1e-3) {w['Ymin']:.1e}, e {w['e']:.1e}{'' if good else '  FAIL'}")
    print("SUCCESS: equilibrium verified" if ok else "FAILURE: equilibrium")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
check_programmed_burn.py: consistency of the programmed-burn speed with the
products equation of state of detonating phases.

    python3 check_programmed_burn.py <case or physicalProperties> ... [--T0 300] [--fix]

For each phase of type `detonating` with JWL or eJWLC products, the CJ point
of the products EOS is computed: the minimum detonation speed on the
products Hugoniot from the reactant state (reactant energy Cv*T0, eConst
reactants with Hf = 0, plus the energy released by the burn: initiation
`E0` [J/m^3] divided by the products rho0, or `e0` [J/kg]). A programmed
burn (`vDet`) slower than this speed has no steady detonation: the energy
released behind the burn front overdrives it and the pressure at the front
grows without bound. A burn faster than it runs on the weak branch (slightly
below p_CJ, stable). With --fix, `vDet` (and `Pcj` if present) of the
initiation dictionary are set to the CJ point of the products EOS.

Products EOS (blastFoam):
    JWL:   p = omega rho e + A (1 - omega/(R1 V)) exp(-R1 V)
                           + B (1 - omega/(R2 V)) exp(-R2 V)
    eJWLC: p = p_s(V) + omega rho (e - e_s(V)),  e_s with the datum E0
with V = rho0/rho.
"""

import argparse
import math
import os
import re
import sys

from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compute_omegaT as co  # noqa: E402


def num(d, key, default=None):
    if key not in d:
        if default is None:
            raise KeyError(key)
        return default
    return co.scalar(d[key])


def cj_point(products, e_r, p0=101325.0):
    eos = products["equationOfState"]
    kind = products["thermoType"]["equationOfState"][0]
    rho0 = num(eos, "rho0")
    v0 = 1/rho0
    A, B, R1, R2, w = (num(eos, k) for k in ("A", "B", "R1", "R2", "omega"))

    if kind == "JWL":
        def p_of(V, e):
            return (w*rho0/V*e + A*(1 - w/(R1*V))*math.exp(-R1*V)
                    + B*(1 - w/(R2*V))*math.exp(-R2*V))
    elif kind == "eJWLC":
        C, E0 = num(eos, "C"), num(eos, "E0")

        def p_of(V, e):
            ps = A*math.exp(-R1*V) + B*math.exp(-R2*V) + C*V**(-(w + 1))
            es = v0*(A/R1*math.exp(-R1*V) + B/R2*math.exp(-R2*V) + C/w*V**(-w)) + E0
            return ps + w*rho0/V*(e - es)
    else:
        return None

    def p_hug(V):
        # e = e_r + (p + p0)/2 (v0 - v); p is linear in e
        v = V*v0
        a = p_of(V, 0.0)
        b = (p_of(V, 1.0) - a)
        return (a + b*(e_r + 0.5*p0*(v0 - v)))/(1 - b*0.5*(v0 - v))

    def speed(V):
        return v0*math.sqrt(max(p_hug(V) - p0, 0.0)/(v0 - V*v0))

    r = minimize_scalar(speed, bounds=(0.3, 0.995), method="bounded", options={"xatol": 1e-10})
    return speed(r.x), p_hug(r.x), rho0/r.x


def detonating_phases(d, path=""):
    for k, v in d.items():
        if isinstance(v, dict):
            if v.get("type") == ["detonating"]:
                yield (path + k), v
            else:
                yield from detonating_phases(v, path + k + "/")


def physical_properties_files(target):
    if os.path.isfile(target):
        return [target]
    out = []
    for root, _, files in os.walk(target):
        for f in files:
            if f.startswith("physicalProperties") and "/processor" not in root:
                out.append(os.path.join(root, f))
    return sorted(out)


def set_entry(text, phase_start, key, value):
    """Set `key value;` in the initiation dictionary following phase_start."""
    i = text.index("initiation", phase_start)
    j = text.index("}", i)
    block = text[i:j]
    new, n = re.subn(rf"(\n\s*{key}\s+)[^;]+;([^\n]*)", rf"\g<1>{value};", block, count=1)
    return text[:i] + new + text[j:], n


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("targets", nargs="+")
    ap.add_argument("--T0", type=float, default=300.0, help="initial reactant temperature [K]")
    ap.add_argument("--fix", action="store_true", help="set vDet (and Pcj) to the CJ point")
    ap.add_argument("--tol", type=float, default=0.005, help="relative tolerance on vDet")
    a = ap.parse_args()

    bad = False
    for target in a.targets:
        for path in physical_properties_files(target):
            try:
                d = co.parse_file(path)
            except Exception:
                continue
            for name, ph in detonating_phases(d):
                ini = ph.get("initiation", {})
                if "vDet" not in ini:
                    continue
                cv = num(ph["reactants"]["thermodynamics"], "Cv")
                rho0p = num(ph["products"]["equationOfState"], "rho0")
                release = num(ini, "E0", 0.0)/rho0p if "E0" in ini else num(ini, "e0", 0.0)
                res = cj_point(ph["products"], cv*a.T0 + release)
                if res is None:
                    continue
                D, p, rho = res
                vdet = num(ini, "vDet")
                err = vdet/D - 1
                flag = ("ok" if abs(err) <= a.tol else
                        "OVERDRIVEN (burn slower than the products CJ speed)" if err < 0 else
                        "weak branch (burn faster than the products CJ speed)")
                print(f"{path} [{name}]: vDet {vdet:.2f} m/s, products CJ {D:.2f} m/s, "
                      f"{p/1e9:.2f} GPa, {rho:.1f} kg/m3 ({100*err:+.2f}%): {flag}")
                if err < -a.tol:
                    bad = True
                if a.fix and abs(err) > a.tol:
                    text = open(path).read()
                    start = text.index(name.split("/")[-1])
                    text, _ = set_entry(text, start, "vDet", f"{D:.1f}")
                    text, _ = set_entry(text, start, "Pcj", f"{p:.4g}")
                    open(path, "w").write(text)
                    print(f"    set vDet {D:.1f}" + (f", Pcj {p:.4g}" if "Pcj" in ini else ""))
    sys.exit(1 if bad and not a.fix else 0)


if __name__ == "__main__":
    main()

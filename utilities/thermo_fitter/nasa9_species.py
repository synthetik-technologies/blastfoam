#!/usr/bin/env python3
"""
nasa9_species.py: species dictionaries with NASA-9 thermodynamics ('thermo
nasa9') for the blastFoam multicomponent gas phases.

All species share the temperature ranges 200 - 1000 - 6000 - 20000 K, as the
nasa9 mixture rule requires. Sources, in order of preference:

  1. EuDetCode gas database (NASA-9, McBride, Zehe and Gordon 2002), e.g.
     external/eudetcode/Data/GasFlame.json: copied exactly when the species
     uses the common ranges.
  2. Cantera NASA-7 data (nasa_gas.yaml, nasa_condensed.yaml, and for
     species found nowhere else gri30.yaml, e.g. CH2O; McBride,
     Gordon and Reno 1993): written exactly in NASA-9 form when the ranges are
     200 - 1000 - 6000 K, with a constant-Cp range 6000 - 20000 K.
  3. Otherwise (other ranges, condensed species): least-squares fit of the
     NASA-9 Cp polynomial in each common range to the source Cp (held
     constant outside the source range), with b1, b2 chosen so that H and S
     are continuous at the range limits and equal to the source at the first
     point. The fit error is reported.

molWeight, the Abel-Nobel covolume b (Tc, Pc) and the Sutherland transport
coefficients are taken from an existing species file (--base); for species
not in it they are copied from a parent species (per mass), see PARENT.

Usage:
    nasa9_species.py --base thermo.compressibleGas --species O2 N2 O N NO ...
        [--eudet external/eudetcode/Data/GasFlame.json]
        [--cantera-gas nasa_gas.yaml] [--cantera-condensed nasa_condensed.yaml]
        [-o thermo.nasa9Gas] [--check]
"""

import argparse
import json
import math
import os
import re
import sys

import numpy as np

T_RANGES = [200.0, 1000.0, 6000.0, 20000.0]
T_STD = 298.15

ATOMIC_WEIGHT = {"H": 1.00794, "C": 12.0107, "N": 14.0067, "O": 15.9994, "Al": 26.981539}

# Cantera names of species not in the EuDetCode gas database
CANTERA_NAME = {"Al": "AL", "AlO": "ALO", "C(gr)": "C(gr)", "Al2O3(s)": "AL2O3(a)",
                "Al2O3(L)": "AL2O3(L)", "Al2O": "AL2O", "AlO2": "ALO2",
                "Al2O2": "AL2O2"}

# Species without an entry in the base file take b, Tc, Pc and transport
# from a parent (b is per mass, so atoms take the value of their molecule)
PARENT = {"O": "O2", "N": "N2", "H": "H2", "OH": "H2O", "NO": "N2", "C": "CO",
          "AlO": "Al", "Al2O": "Al", "AlO2": "Al", "Al2O2": "Al",
          "Al2O3(L)": "Al2O3(s)", "CH": "CH4", "NH": "NH3", "HO2": "O2", "H2O2": "H2O",
          "CH3": "CH4", "CH2O": "CH4", "HCO": "CO", "CH3O": "CH4"}


# --- NASA polynomials ------------------------------------------------------

def cp9(a, T):
    return a[0]/T**2 + a[1]/T + a[2] + T*(a[3] + T*(a[4] + T*(a[5] + T*a[6])))


def h9(a, T):
    """H/R"""
    return (-a[0]/T + a[1]*math.log(T) + a[7]
            + T*(a[2] + T*(a[3]/2 + T*(a[4]/3 + T*(a[5]/4 + T*a[6]/5)))))


def s9(a, T):
    """S/R"""
    return (-a[0]/(2*T**2) - a[1]/T + a[2]*math.log(T) + a[8]
            + T*(a[3] + T*(a[4]/2 + T*(a[5]/3 + T*a[6]/4))))


def nasa7_to_9(c):
    """NASA-7 (a1..a5, b1, b2) -> NASA-9 (0, 0, a1..a5, b1, b2), exact."""
    return [0.0, 0.0, c[0], c[1], c[2], c[3], c[4], c[5], c[6]]


class Piecewise:
    """Cp/R, H/R, S/R from a list of (Tmin, Tmax, nasa9 coeffs); constant Cp
    outside the overall range."""

    def __init__(self, ranges):
        self.ranges = ranges
        self.Tmin = ranges[0][0]
        self.Tmax = ranges[-1][1]

    def coeffs(self, T):
        # The upper range at a range limit, as nasa9Thermo
        for lo, hi, a in self.ranges:
            if T < hi:
                return a
        return self.ranges[-1][2]

    def cp(self, T):
        Tf = min(max(T, self.Tmin), self.Tmax)
        return cp9(self.coeffs(Tf), Tf)

    def h(self, T):
        Tf = min(max(T, self.Tmin), self.Tmax)
        return h9(self.coeffs(Tf), Tf) + cp9(self.coeffs(Tf), Tf)*(T - Tf)

    def s(self, T):
        Tf = min(max(T, self.Tmin), self.Tmax)
        return s9(self.coeffs(Tf), Tf) + cp9(self.coeffs(Tf), Tf)*math.log(T/Tf)


# --- Sources ---------------------------------------------------------------

def load_eudet(path):
    if not path or not os.path.isfile(path):
        return {}
    return json.load(open(path))["species"]


def load_cantera(paths):
    import yaml
    out = {}
    for p in paths:
        if p and os.path.isfile(p):
            for sp in yaml.safe_load(open(p))["species"]:
                out[sp["name"]] = sp
    return out


def continuous_ranges(Ts, cs):
    """Number of ranges from the start that join continuously (Cp within
    1%, H/R within 5 K, S/R within 0.01); EuDetCode has a discontinuous
    6000 - 20000 K range for some polyatomic species (H2O, NH3, CH4)."""
    for i in range(1, len(cs)):
        T, a, b = Ts[i], cs[i - 1], cs[i]
        if (abs(cp9(a, T) - cp9(b, T)) > 0.01*abs(cp9(a, T))
                or abs(h9(a, T) - h9(b, T)) > 5 or abs(s9(a, T) - s9(b, T)) > 0.01):
            return i
    return len(cs)


def source(name, eudet, cantera):
    """Return (Piecewise, origin)."""
    if name in eudet:
        n = eudet[name]["nasaParams"]
        Ts, cs = n["Tswitch"], n["coeffs"]
        k = continuous_ranges(Ts, cs)
        origin = "EuDetCode NASA-9"
        if k < len(cs):
            origin += f" (range above {Ts[k]:g} K discontinuous, dropped)"
        return Piecewise([(Ts[i], Ts[i + 1], list(cs[i])) for i in range(k)]), origin
    cname = CANTERA_NAME.get(name, name)
    if cname in cantera:
        th = cantera[cname]["thermo"]
        Ts = th["temperature-ranges"]
        data = th["data"]
        if th["model"] == "NASA7":
            return Piecewise([(Ts[i], Ts[i + 1], nasa7_to_9(data[i])) for i in range(len(data))]), f"Cantera NASA-7 ({cname})"
        if th["model"] == "NASA9":
            return Piecewise([(Ts[i], Ts[i + 1], list(data[i])) for i in range(len(data))]), f"Cantera NASA-9 ({cname})"
    sys.exit(f"No thermodynamic data for {name}")


def on_common_ranges(src):
    """NASA-9 coefficients on T_RANGES, and the max relative Cp/H error."""
    ranges = [(lo, hi) for lo, hi, _ in src.ranges]
    common = list(zip(T_RANGES[:-1], T_RANGES[1:]))
    out = []
    if ranges == common:
        return [list(a) for _, _, a in src.ranges], 0.0
    if ranges == common[:2]:
        out = [list(a) for _, _, a in src.ranges]
        # Constant Cp above 6000 K, H and S continuous
        T = T_RANGES[2]
        cp = src.cp(T)
        out.append([0, 0, cp, 0, 0, 0, 0, src.h(T) - cp*T, src.s(T) - cp*math.log(T)])
        return out, 0.0
    # Least-squares fit of Cp/R in each common range (columns scaled, the
    # basis spans T^-2 to T^4); a range beyond the source data has the
    # constant Cp of the source exactly
    for lo, hi in common:
        if lo >= src.Tmax:
            out.append([0.0, 0.0, src.cp(src.Tmax), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            continue
        T = np.geomspace(lo, hi, 400)
        A = np.stack([T**-2, T**-1, np.ones_like(T), T, T**2, T**3, T**4], axis=1)
        scale = np.linalg.norm(A, axis=0)
        y = np.array([src.cp(t) for t in T])
        a = list(np.linalg.lstsq(A/scale, y, rcond=None)[0]/scale)
        out.append(a + [0.0, 0.0])
    # Integration constants: continuous with the source at 200 K, then
    # continuous between the ranges
    T0 = T_RANGES[0]
    out[0][7] = src.h(T0) - h9(out[0], T0)
    out[0][8] = src.s(T0) - s9(out[0], T0)
    for i in range(1, len(out)):
        T = T_RANGES[i]
        out[i][7] = h9(out[i - 1], T) - h9(out[i], T)
        out[i][8] = s9(out[i - 1], T) - s9(out[i], T)
    err = 0.0
    fit = Piecewise([(lo, hi, a) for (lo, hi), a in zip(common, out)])
    for T in np.geomspace(T_RANGES[0], T_RANGES[-1], 400):
        err = max(err, abs(fit.cp(T)/src.cp(T) - 1))
    return out, err


# --- Base species file -----------------------------------------------------

def parse_base(path):
    """Species entries of an OpenFOAM thermo file as raw text blocks."""
    text = open(path).read()
    blocks = {}
    for m in re.finditer(r'^("?[A-Za-z0-9()]+"?)\s*\n\{', text, re.M):
        name = m.group(1).strip('"')
        depth, i = 0, m.end() - 1
        for j in range(i, len(text)):
            if text[j] == "{":
                depth += 1
            elif text[j] == "}":
                depth -= 1
                if depth == 0:
                    blocks[name] = text[i:j + 1]
                    break
    return blocks


def sub(block, name):
    m = re.search(name + r"\s*\{([^}]*)\}", block)
    return m.group(1) if m else ""


def molweight(name):
    toks = re.findall(r"([A-Z][a-z]?)(\d*)", name.split("(")[0])
    return sum(ATOMIC_WEIGHT[e]*(int(n) if n else 1) for e, n in toks)


# --- Output ----------------------------------------------------------------

def entry(name, coeffs, base, origin, err):
    blk = base.get(name) or base.get(PARENT.get(name, ""), "")
    if not blk:
        sys.exit(f"No base entry for {name} or its parent {PARENT.get(name)}")
    W = molweight(name)
    eos = "\n".join("        " + l.strip() for l in sub(blk, "equationOfState").strip().splitlines())
    tr = "\n".join("        " + l.strip() for l in sub(blk, "transport").strip().splitlines())
    rows = "\n".join("            (" + " ".join(f"{c:.10g}" for c in a) + ")" for a in coeffs)
    note = f"// {origin}" + (f"; fitted, max Cp error {100*err:.2f}%" if err > 0 else "")
    if name not in base:
        note += f"; b and transport from {PARENT[name]}"
    q = f'"{name}"' if "(" in name else name
    return f"""{q}
{{
    {note}
    specie
    {{
        molWeight       {W:.10g};
    }}
    equationOfState
    {{
{eos}
    }}
    thermodynamics
    {{
        Tswitch         ({" ".join(f"{t:g}" for t in T_RANGES)});
        coeffs
        (
{rows}
        );
    }}
    transport
    {{
{tr}
    }}
}}
"""


def main():
    root = "/opt/OpenFOAM/blastfoam"
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--base", required=True)
    ap.add_argument("--species", nargs="+", required=True)
    ap.add_argument("--eudet", default=os.path.join(root, "external/eudetcode/Data/GasFlame.json"))
    ap.add_argument("--cantera-gas", default="/usr/local/share/cantera/data/nasa_gas.yaml")
    ap.add_argument("--cantera-mech", default="/usr/local/share/cantera/data/gri30.yaml",
                    help="lowest-priority source (GRI-3.0 NASA-7, e.g. for CH2O)")
    ap.add_argument("--cantera-condensed", default="/usr/local/share/cantera/data/nasa_condensed.yaml")
    ap.add_argument("-o", "--output")
    ap.add_argument("--check", action="store_true", help="compare Cp, H, S with the source")
    args = ap.parse_args()

    eudet = load_eudet(args.eudet)
    cantera = load_cantera([args.cantera_mech, args.cantera_gas, args.cantera_condensed])
    base = parse_base(args.base)

    out = ["// Generated by utilities/thermo_fitter/nasa9_species.py\n"]
    for name in args.species:
        src, origin = source(name, eudet, cantera)
        coeffs, err = on_common_ranges(src)
        out.append(entry(name, coeffs, base, origin, err))
        if args.check:
            fit = Piecewise([(lo, hi, a) for lo, hi, a in zip(T_RANGES[:-1], T_RANGES[1:], coeffs)])
            e = [max(abs(fit.cp(T) - src.cp(T))/abs(src.cp(T)),
                     abs(fit.h(T) - src.h(T))/max(abs(src.h(T)), 1e3),
                     abs(fit.s(T) - src.s(T))/abs(src.s(T)))
                 for T in np.geomspace(T_RANGES[0], T_RANGES[-1], 500)]
            print(f"{name:10s} {origin:32s} max rel. diff Cp/H/S vs source: {max(e):.2e}", file=sys.stderr)
    text = "\n".join(out)
    if args.output:
        open(args.output, "w").write(text)
    else:
        print(text)


if __name__ == "__main__":
    main()

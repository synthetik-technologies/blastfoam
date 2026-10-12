#!/usr/bin/env python3
"""
compute_omegaT.py: Compute the eJWLC isentrope-temperature exponent omegaT
for a blastFoam case using pressurePhaseChange.

The eJWLC isentrope temperature is T_s(V) = T_CJ (V_CJ/V)^omegaT. omegaT is
chosen so that T_s passes through both the CJ temperature (EuDetCode) and the
gas-EOS temperature at the phase-transfer state:

    V_trans : p_s(V_trans) = p_trans  (products isentrope)
    T_trans = p_trans (v_trans - b_mix)/R_mix   (Abel-Nobel gas)
            = p_trans v_trans/R_mix             (perfect/ideal gas)
    omegaT  = ln(T_CJ/T_trans)/ln(V_trans/V_CJ)

with v_trans = V_trans/rho0 and the gas mixture (R_mix, mass-fraction
weighted b_mix) taken from the transfer targetFractions in fvModels.
omegaT changes only the energy-temperature relation, not the pressure
(report section "Energy Conservation and Temperature Consistency During
Transfer").

Usage:
    python3 compute_omegaT.py <caseDir> [--ptrans PA] [--gas-case DIR] [--write]

--inert Y:RHO declares that the eJWLC fit includes a mass fraction Y of
inert condensed material of density RHO (e.g. 0.27:2700 for a Comp-B/Al
fit with inert Al); the gas specific volume is then
v_gas = (v - Y/RHO)/(1 - Y).

--gas-case takes the transfer composition and gas equation of state from
another case, for cases whose own transfer target is a placeholder or a
different gas (e.g. disabled transfer, air as the gas phase).
"""

import argparse
import math
import os
import re
import sys

R_UNIVERSAL = 8314.462618  # J/(kmol K)


# --------------------------------------------------------------------------
# Minimal OpenFOAM dictionary parser (dicts, scalars, words, lists, #include)
# --------------------------------------------------------------------------

def _tokenize(text):
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    text = re.sub(r"//[^\n]*", " ", text)
    # Words may contain parentheses, e.g. specie names C(gr), Al2O3(s)
    word = r"[^\s{}();]+(?:\([^\s{}();]*\)[^\s{}();]*)*"
    return re.findall(r'"[^"]*"|[{};]|' + word + r'|[()]', text)


def parse_file(path):
    with open(path) as f:
        tokens = _tokenize(f.read())
    d, _ = _parse_dict(tokens, 0, os.path.dirname(path))
    return d


def _parse_dict(tokens, i, base):
    d = {}
    while i < len(tokens):
        tok = tokens[i]
        if tok == "}":
            return d, i + 1
        if tok == "#include":
            inc = os.path.join(base, tokens[i + 1].strip('"'))
            d.update(parse_file(inc))
            i += 2
            continue
        if tok.startswith("#"):
            i += 2
            continue
        key = tok.strip('"')
        i += 1
        if i < len(tokens) and tokens[i] == "{":
            sub, i = _parse_dict(tokens, i + 1, base)
            d[key] = sub
            continue
        vals = []
        depth = 0
        while i < len(tokens) and not (tokens[i] == ";" and depth == 0):
            if tokens[i] == "(":
                depth += 1
            elif tokens[i] == ")":
                depth -= 1
            vals.append(tokens[i])
            i += 1
        d[key] = vals
        i += 1
    return d, i


def scalar(entry):
    return float(entry[0]) if isinstance(entry, list) else float(entry)


# --------------------------------------------------------------------------
# Case inspection
# --------------------------------------------------------------------------

def phase_dict(constant, phase):
    """Phase block: physicalProperties.<phase> or <phase> {} in
    physicalProperties."""
    sep = os.path.join(constant, "physicalProperties." + phase)
    if os.path.isfile(sep):
        return parse_file(sep)
    return parse_file(os.path.join(constant, "physicalProperties"))[phase]


def find_ejwlc(d):
    """Return the dictionary holding thermoType{equationOfState eJWLC}."""
    tt = d.get("thermoType")
    if isinstance(tt, dict) and tt.get("equationOfState") == ["eJWLC"]:
        return d
    for v in d.values():
        if isinstance(v, dict):
            r = find_ejwlc(v)
            if r is not None:
                return r
    return None


def transfer_model(constant):
    for name in ("fvModels", "fvModels.disabled"):
        path = os.path.join(constant, name)
        if os.path.isfile(path):
            for v in parse_file(path).values():
                if isinstance(v, dict) and v.get("type", [""])[0] in (
                    "pressurePhaseChange", "pressurePhaseTransfer"
                ):
                    return v, path
    sys.exit("No pressurePhaseChange model found in constant/fvModels*")


def gas_target(model, constant):
    """Name and mass fractions of the gas target (the target whose thermo
    is not eJWLC/condensed; the one with targetFractions if several)."""
    if "targets" in model:
        for name, t in model["targets"].items():
            if "targetFractions" in t:
                return name, t["targetFractions"]
        return next(iter(model["targets"])), {}
    return model["targetPhase"][0], model.get("targetFractions", {})


def p_trans_from(model, default):
    if "pInterval" in model:
        lo, hi = (float(x) for x in model["pInterval"] if x not in "()")
        return 0.5*(lo + hi)
    for lo, hi in (("pMin", "pMax"), ("pLower", "pUpper")):
        if lo in model and hi in model:
            return 0.5*(scalar(model[lo]) + scalar(model[hi]))
    if "pThreshold" in model:
        return scalar(model["pThreshold"])
    return default


def gas_mixture(gas, fractions):
    """Return (R_mix, b_mix, eos) for the gas phase at the given mass
    fractions (b mass-fraction weighted, as in AbelNobel)."""
    eos = gas["thermoType"]["equationOfState"][0]
    if not fractions:
        # Single-specie gas defined inline
        W = scalar(gas["specie"]["molWeight"])
        b = scalar(gas.get("equationOfState", {}).get("b", [0.0]))
        return R_UNIVERSAL/W, b, eos
    Ysum = sum(scalar(v) for v in fractions.values())
    R_mix = b_mix = 0.0
    for sp, y in fractions.items():
        y = scalar(y)/Ysum
        data = gas[sp]
        R_mix += y*R_UNIVERSAL/scalar(data["specie"]["molWeight"])
        b_mix += y*scalar(data.get("equationOfState", {}).get("b", [0.0]))
    return R_mix, b_mix, eos


def isentrope_volume(eos, p):
    """Solve p_s(V) = p for the relative volume V (bisection)."""
    A, B, C = (scalar(eos[k]) for k in ("A", "B", "C"))
    R1, R2, w = (scalar(eos[k]) for k in ("R1", "R2", "omega"))

    def ps(V):
        return A*math.exp(-R1*V) + B*math.exp(-R2*V) + C*V**(-(w + 1.0))

    lo, hi = 1e-3, 1e5
    for _ in range(200):
        mid = math.sqrt(lo*hi)
        if ps(mid) > p:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo*hi)


def write_omegaT(path, value):
    """Insert or update omegaT in the eJWLC equationOfState block."""
    with open(path) as f:
        text = f.read()
    m = re.search(r"equationOfState\s+eJWLC;", text)
    blk = re.search(r"equationOfState\s*\{", text[m.end():])
    start = m.end() + blk.end()
    end = text.index("}", start)
    body = text[start:end]
    line = f"omegaT          {value:.6f};"
    if re.search(r"\bomegaT\b", body):
        body = re.sub(r"omegaT\s+[-0-9.eE+]+;", line, body)
    else:
        anchor = re.search(r"\n(\s*)(E0|rho_CJ)\s+[^\n]*;[^\n]*\n", body)
        if anchor is None:
            sys.exit(f"Cannot find E0/rho_CJ in eJWLC block of {path}")
        indent = anchor.group(1)
        body = (
            body[:anchor.end()]
          + f"{indent}{line} // T_s through T_CJ and the gas-EOS T at transfer"
          + "\n" + body[anchor.end():]
        )
    with open(path, "w") as f:
        f.write(text[:start] + body + text[end:])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("case")
    ap.add_argument("--ptrans", type=float, default=1.0e7,
                    help="Transfer pressure [Pa] if not given in fvModels")
    ap.add_argument("--gas-case", default=None,
                    help="Case providing the transfer composition and gas "
                         "EOS (default: the case itself)")
    ap.add_argument("--inert", default=None, metavar="Y:RHO",
                    help="Inert mass fraction and density included in the "
                         "eJWLC fit, e.g. 0.27:2700")
    ap.add_argument("--write", action="store_true",
                    help="Write omegaT into the products eJWLC block")
    args = ap.parse_args()

    constant = os.path.join(args.case, "constant")
    model, model_path = transfer_model(constant)
    source = model["sourcePhase"][0]
    p_trans = p_trans_from(model, args.ptrans)

    gas_constant = constant
    if args.gas_case:
        gas_constant = os.path.join(args.gas_case, "constant")
        model_g, _ = transfer_model(gas_constant)
        gas_name, fractions = gas_target(model_g, gas_constant)
    else:
        gas_name, fractions = gas_target(model, constant)

    src = phase_dict(constant, source)
    prod = find_ejwlc(src)
    if prod is None:
        sys.exit(f"Source phase {source} has no eJWLC equation of state")
    eos = prod["equationOfState"]
    rho0 = scalar(eos["rho0"])
    T_CJ = scalar(eos["T_CJ"])
    V_CJ = rho0/scalar(eos["rho_CJ"])
    omega = scalar(eos["omega"])

    V_trans = isentrope_volume(eos, p_trans)
    v_trans = V_trans/rho0
    if args.inert:
        Y_inert, rho_inert = (float(x) for x in args.inert.split(":"))
        v_trans = (v_trans - Y_inert/rho_inert)/(1.0 - Y_inert)

    R_mix, b_mix, gas_eos = gas_mixture(
        phase_dict(gas_constant, gas_name), fractions
    )
    if gas_eos.startswith("AbelNobel"):
        T_trans = p_trans*(v_trans - b_mix)/R_mix
    elif gas_eos in ("perfectGas", "idealGas"):
        T_trans = p_trans*v_trans/R_mix
    else:
        sys.exit(f"Unsupported gas equation of state {gas_eos}")

    omegaT = math.log(T_CJ/T_trans)/math.log(V_trans/V_CJ)
    Ts_omega = T_CJ*(V_CJ/V_trans)**omega

    print(f"case            {args.case}")
    print(f"transfer model  {os.path.basename(model_path)}: "
          f"{source} -> {gas_name} ({gas_eos}) at {p_trans/1e5:.1f} bar")
    if args.gas_case:
        print(f"gas from        {args.gas_case}")
    print(f"V_trans         {V_trans:.4f} (rho {rho0/V_trans:.4f} kg/m3)")
    if args.inert:
        print(f"inert           Y = {Y_inert}, rho = {rho_inert} kg/m3; "
              f"gas v = {v_trans:.5e} m3/kg")
    print(f"R_mix, b_mix    {R_mix:.3f} J/(kg K), {b_mix:.4e} m3/kg")
    print(f"T_trans         {T_trans:.2f} K "
          f"(power law with omega gives {Ts_omega:.2f} K)")
    print(f"omega, omegaT   {omega:.4f}, {omegaT:.6f}")

    if args.write:
        files = [
            os.path.join(constant, "physicalProperties." + source),
            os.path.join(constant, "physicalProperties"),
        ]
        path = next(f for f in files if os.path.isfile(f))
        write_omegaT(path, omegaT)
        print(f"wrote omegaT to {path}")


if __name__ == "__main__":
    main()

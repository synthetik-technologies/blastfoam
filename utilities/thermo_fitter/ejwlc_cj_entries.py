#!/usr/bin/env python3
"""
ejwlc_cj_entries.py: eJWLC products entries and programmed-burn initiation of
a detonating phase from an EuDetCode CJDetonation run, without phase transfer
(e.g. tests/1D/planar_znd).

    python3 ejwlc_cj_entries.py <eudet_dir> [--e-reactants 309000]

<eudet_dir> holds cj_output.txt, jwl_output.txt and isentrope_output.dat of
    CJDetonation "<explosive>" <rho0> BKW -o <eudet_dir>

The JWL fit of EuDetCode is used as the eJWLC reference isentrope
p_s(V) = A exp(-R1 V) + B exp(-R2 V) + C V^-(omega+1). Its energy datum E0 is
set so that the isentrope passes through the CJ point on the products
Hugoniot of the reactants as the case defines them,
    e_s(V_CJ) = e_reactants + (p_CJ + p0)/2 (v0 - v_CJ),
with e_reactants the reactant internal energy at the initial state in the
case's datum (Cv*T0 for an eConst reactant with Hf = 0). The CJ point of
the resulting products EOS (minimum detonation speed on its Hugoniot) is
then checked against D_CJ: a products EOS whose own CJ speed differs from
the programmed-burn speed has no steady detonation at that speed (an energy
release that is too large overdrives the front and the pressure behind the
burn front grows without bound). omegaT is fitted to the EuDetCode isentrope
temperatures, T_s = T_CJ (V_CJ/V)^omegaT, between the CJ point and 1 GPa.
"""

import argparse
import math
import os
import sys

import numpy as np
from scipy.optimize import minimize_scalar


def parse(path):
    d = {}
    for line in open(path):
        if ":" in line:
            k, v = line.split(":", 1)
            try:
                d[k.strip()] = float(v.split()[0])
            except (ValueError, IndexError):
                pass
    return d


def compute(eudet_dir, e_reactants=1030.0*300.0, p0=101325.0):
    cj = parse(os.path.join(eudet_dir, "cj_output.txt"))
    jwl = parse(os.path.join(eudet_dir, "jwl_output.txt"))
    rho0 = cj["Initial Density"] if "Initial Density" in cj else None
    if rho0 is None:
        for line in open(os.path.join(eudet_dir, "cj_output.txt")):
            if line.startswith("Initial Density"):
                rho0 = float(line.split(":")[1].split()[0])
    p_CJ = cj["CJ Pressure (Pa)"]
    T_CJ = cj["CJ Temperature (K)"]
    v_CJ = cj["CJ Volume (m^3/kg)"]
    D = cj["CJ Velocity (m/s)"]
    A, B, C = (jwl[k]*1e9 for k in ("A", "B", "C"))
    R1, R2, w = jwl["R1"], jwl["R2"], jwl["omega"]
    v0 = 1/rho0
    V_CJ = v_CJ/v0

    def ps(V):
        return A*math.exp(-R1*V) + B*math.exp(-R2*V) + C*V**(-(w + 1))

    def eb(V):
        return v0*(A/R1*math.exp(-R1*V) + B/R2*math.exp(-R2*V) + C/w*V**(-w))

    E0 = e_reactants + 0.5*(p_CJ + p0)*(v0 - v_CJ) - eb(V_CJ)

    # CJ point of the products EOS p = ps + w rho (e - es) on the Hugoniot
    def p_hug(V):
        v = V*v0
        rho = 1/v
        es = eb(V) + E0
        return (ps(V) + w*rho*(e_reactants + 0.5*p0*(v0 - v) - es))/(1 - w*rho*0.5*(v0 - v))

    def speed(V):
        return v0*math.sqrt(max(p_hug(V) - p0, 0)/(v0 - V*v0))

    r = minimize_scalar(speed, bounds=(0.3, 0.99), method="bounded", options={"xatol": 1e-10})
    D_eos, p_eos, V_eos = speed(r.x), p_hug(r.x), r.x

    # omegaT from the isentrope temperatures between CJ and 1 GPa
    data = np.loadtxt(os.path.join(eudet_dir, "isentrope_output.dat"))
    sel = data[:, 1] >= 1e9
    slope = np.polyfit(np.log(data[sel, 0]/v_CJ), np.log(data[sel, 2]/T_CJ), 1)[0]
    omegaT = -slope

    return dict(rho0=rho0, A=A, B=B, C=C, R1=R1, R2=R2, omega=w, Cv=jwl["Cv"],
                T_CJ=T_CJ, rho_CJ=1/v_CJ, p_CJ=p_CJ, D=D, E0=E0, omegaT=omegaT,
                D_eos=D_eos, p_eos=p_eos, rho_eos=rho0/V_eos)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("eudet_dir")
    ap.add_argument("--e-reactants", type=float, default=1030.0*300.0,
                    help="reactant internal energy at the initial state [J/kg]")
    a = ap.parse_args()
    r = compute(a.eudet_dir, a.e_reactants)
    print(f"// EuDetCode CJ: D {r['D']:.2f} m/s, p_CJ {r['p_CJ']/1e9:.4f} GPa, rho_CJ {r['rho_CJ']:.2f}, T_CJ {r['T_CJ']:.2f} K")
    print(f"// CJ point of the products EOS: D {r['D_eos']:.2f} m/s ({100*(r['D_eos']/r['D'] - 1):+.3f}%), "
          f"p {r['p_eos']/1e9:.4f} GPa ({100*(r['p_eos']/r['p_CJ'] - 1):+.2f}%), rho {r['rho_eos']:.1f}")
    print("products equationOfState")
    print("{")
    for k, fmt in (("rho0", ".2f"), ("A", ".6g"), ("B", ".6g"), ("C", ".6g"), ("R1", ".6g"),
                   ("R2", ".6g"), ("omega", ".6g"), ("Cv", ".6g"), ("T_CJ", ".2f"),
                   ("rho_CJ", ".2f"), ("omegaT", ".6f"), ("E0", ".2f")):
        print(f"    {k:15s} {r[k]:{fmt}};")
    print("}")
    print(f"initiation: vDet {r['D']:.2f}; Pcj {r['p_CJ']:.6g};")
    if abs(r["D_eos"]/r["D"] - 1) > 5e-3:
        sys.exit("The products EOS is not consistent with D_CJ (> 0.5%)")


if __name__ == "__main__":
    main()

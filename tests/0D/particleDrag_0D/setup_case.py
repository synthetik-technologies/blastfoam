#!/usr/bin/env python3
"""Initial fields of the 0D drag case: gas at rest, particles moving.

    python3 setup_case.py --d 1e-4 --p 1e5 --T 300 --U 600 --alpha 1e-6 \
        --deltaT 1e-7 --endTime 2e-3
"""

import argparse
import os

CASE = os.path.dirname(os.path.abspath(__file__))

a = argparse.ArgumentParser()
for k in ("d", "p", "T", "U", "alpha", "deltaT", "endTime"):
    a.add_argument("--" + k, type=float, required=True)
a = a.parse_args()

# Ideal gas (eConst, Cv 718, gamma 1.4) and Murnaghan Al as in constant/
rho_g = a.p/(0.4*718.0*a.T)
rho_al = 2700.0*(1.0 + 4.3*(a.p - 101325.0)/72.7e9)**(1.0/4.3)

template = open(os.path.join(CASE, "0.orig", "T")).read()


def write(name, cls, dims, val):
    s = (template.replace("class       volScalarField;", f"class       {cls};")
         .replace("object      T;", f"object      {name};")
         .replace("dimensions      [0 0 0 1 0 0 0];", f"dimensions      {dims};")
         .replace("internalField   uniform TVALUE;", f"internalField   uniform {val};"))
    open(os.path.join(CASE, "0", name), "w").write(s)


os.makedirs(os.path.join(CASE, "0"), exist_ok=True)
s, v, t = "volScalarField", "volVectorField", "[0 0 0 1 0 0 0]"
write("alpha.Al", s, "[0 0 0 0 0 0 0]", a.alpha)
write("alpha.gas", s, "[0 0 0 0 0 0 0]", 1 - a.alpha)
write("p", s, "[1 -1 -2 0 0 0 0]", a.p)
for n in ("T", "T.gas", "T.Al"):
    write(n, s, t, a.T)
write("rho.gas", s, "[1 -3 0 0 0 0 0]", f"{rho_g:.10g}")
write("rho.Al", s, "[1 -3 0 0 0 0 0]", f"{rho_al:.10g}")
write("Theta.Al", s, "[0 2 -2 0 0 0 0]", 0)
write("U", v, "[0 1 -1 0 0 0 0]", "(0 0 0)")
write("U.gas", v, "[0 1 -1 0 0 0 0]", "(0 0 0)")
write("U.Al", v, "[0 1 -1 0 0 0 0]", f"({a.U} 0 0)")

for f, key, val in (("constant/phaseProperties", "D_PARTICLE", a.d),
                    ("system/controlDict", "DELTA_T", a.deltaT),
                    ("system/controlDict", "END_TIME", a.endTime)):
    path = os.path.join(CASE, f)
    text = open(path).read()
    open(path, "w").write(text.replace(key, repr(val)))

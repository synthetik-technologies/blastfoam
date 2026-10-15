#!/usr/bin/env python3
"""Initial fields of the 0D Al burning case: particles and air at rest at
the same temperature (above the ignition temperature).

    python3 setup_case.py --d 1e-5 --p 101325 --T 2500 --alpha 1e-7 \
        --deltaT 1e-7 --endTime 1e-3
"""

import argparse
import os

CASE = os.path.dirname(os.path.abspath(__file__))

a = argparse.ArgumentParser()
for k in ("d", "p", "T", "alpha", "deltaT", "endTime"):
    a.add_argument("--" + k, type=float, required=True)
a = a.parse_args()

# Air (mass fractions) as an Abel-Nobel gas with the co-volumes of
# constant/thermo.nasa9Gas, and Murnaghan Al as in constant/
Y = {"O2": 0.2331, "N2": 0.7669, "Al": 0.0}
W = {"O2": 31.9988, "N2": 28.0134, "Al": 26.9815}
thermo = open(os.path.join(CASE, "constant", "thermo.nasa9Gas")).read()
b = {}
for sp in Y:
    blk = thermo[thermo.index("\n" + sp + "\n"):]
    b[sp] = float(blk.split("b ", 1)[1].split(";")[0])
R = 8314.462618/(1.0/sum(Y[k]/W[k] for k in Y))
bmix = sum(Y[k]*b[k] for k in Y)
rho_g = 1.0/(R*a.T/a.p + bmix)
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
write("d.Al", s, "[0 1 0 0 0 0 0]", a.d)
write("Theta.Al", s, "[0 2 -2 0 0 0 0]", 0)
write("U", v, "[0 1 -1 0 0 0 0]", "(0 0 0)")
write("U.gas", v, "[0 1 -1 0 0 0 0]", "(0 0 0)")
write("U.Al", v, "[0 1 -1 0 0 0 0]", "(0 0 0)")
for sp, y in Y.items():
    write(sp + ".gas", s, "[0 0 0 0 0 0 0]", y)
write("Ydefault.gas", s, "[0 0 0 0 0 0 0]", 0)

for f, key, val in (("constant/phaseProperties", "D_PARTICLE", a.d),
                    ("system/controlDict", "DELTA_T", a.deltaT),
                    ("system/controlDict", "END_TIME", a.endTime)):
    path = os.path.join(CASE, f)
    text = open(path).read()
    open(path, "w").write(text.replace(key, repr(val)))

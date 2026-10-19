#!/usr/bin/env python3
"""Initial fields of the 0D Al burning case with products: particles and
gas at rest at the same temperature.

    python3 setup_case.py --d 1e-5 --p 101325 --T 2500 --alpha 1e-7 \
        --deltaT 1e-7 --endTime 1e-3 --Y O2=0.2331,N2=0.7669 \
        --Tign 2000 --Textinction 0 [--Tp 1500] [--Kh 0] [--ignitedFraction false]
"""

import argparse
import os
import re

CASE = os.path.dirname(os.path.abspath(__file__))

a = argparse.ArgumentParser()
for k in ("d", "p", "T", "alpha", "deltaT", "endTime", "Tign", "Textinction"):
    a.add_argument("--" + k, type=float, required=True)
a.add_argument("--Y", required=True)
a.add_argument("--Tp", type=float, default=None, help="particle temperature (default T)")
a.add_argument("--Kh", type=float, default=0.0, help="heat transfer coefficient [W/m3/K]")
a.add_argument("--ignitedFraction", default="false")
a = a.parse_args()

# Air (mass fractions) as an Abel-Nobel gas with the co-volumes of
# constant/thermo.nasa9Gas, and Murnaghan Al as in constant/
SPECIES = ["O2", "N2", "H2O", "CO2", "H2", "CO", "Al2O3(s)", "Al2O3(L)"]
Y = dict.fromkeys(SPECIES, 0.0)
for kv in a.Y.split(","):
    k, v = kv.split("=")
    Y[k] = float(v)
W = {"O2": 31.9988, "N2": 28.0134, "H2O": 18.01528, "CO2": 44.0095, "H2": 2.01588,
     "CO": 28.0101, "Al2O3(s)": 101.96128, "Al2O3(L)": 101.96128}
thermo = open(os.path.join(CASE, "constant", "thermo.nasa9Gas")).read()
b = {}
for sp in Y:
    key = '"' + sp + '"' if "(" in sp else sp
    blk = thermo[thermo.index("\n" + key + "\n"):]
    b[sp] = float(re.search(r"\n\s+b\s+([-0-9.eE+]+);", blk).group(1))
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
for n in ("T", "T.gas"):
    write(n, s, t, a.T)
write("T.Al", s, t, a.Tp if a.Tp is not None else a.T)
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
                    ("system/controlDict", "END_TIME", a.endTime),
                    ("constant/phaseProperties", "TIGN", a.Tign),
                    ("constant/phaseProperties", "TEXTINCTION", a.Textinction),
                    ("constant/phaseProperties", "KH", a.Kh),
                    ("constant/phaseProperties", "IGNITEDFRACTION", a.ignitedFraction)):
    path = os.path.join(CASE, f)
    text = open(path).read()
    open(path, "w").write(text.replace(key, val if isinstance(val, str) else repr(val)))

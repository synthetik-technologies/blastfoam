#!/usr/bin/env python3
"""
setup_case.py: initial fields and time step of the 0D gas kinetics case
(uniform T, p, zero velocity, given mass fractions).

    python3 setup_case.py --T 1500 --p 1e6 --Y "CO:0.14,O2:0.12,N2:0.74" \
        --deltaT 1e-8 --endTime 1e-4
"""

import argparse
import os

SPECIES = ["O2", "N2", "CO", "H2O", "CO2", "NH3", "CH4", "C2H6", "Al", "H2",
           "O", "N", "NO", "H", "OH", "AlO", "HO2", "H2O2", "CH3", "CH2O", "HCO",
           "CH3O", "C(gr)", "Al2O3(s)"]
def mixture(spec):
    Y = {}
    for item in spec.split(","):
        k, v = item.split(":")
        Y[k.strip()] = float(v)
    unknown = set(Y) - set(SPECIES)
    if unknown:
        raise SystemExit(f"Unknown species {sorted(unknown)}")
    s = sum(Y.values())
    return {k: Y.get(k, 0)/s for k in SPECIES}


HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  12
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       {cls};
    location    "0";
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      {dim};

internalField   uniform {val};

boundaryField
{{
    walls
    {{
        type            {wall};
    }}
    defaultFaces
    {{
        type            empty;
    }}
}}
"""


def write(name, val, dim, cls="volScalarField", wall="zeroGradient"):
    obj = f'"{name}"' if "(" in name else name
    with open(os.path.join("0", name), "w") as f:
        f.write(HEADER.format(cls=cls, obj=obj, dim=dim, val=val, wall=wall))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=float, default=4000)
    ap.add_argument("--p", type=float, default=1e7)
    ap.add_argument("--Y", required=True, help='mass fractions "Sp:Y,..."')
    ap.add_argument("--deltaT", type=float, required=True)
    ap.add_argument("--endTime", type=float, required=True)
    a = ap.parse_args()

    os.makedirs("0", exist_ok=True)
    write("T", repr(a.T), "[0 0 0 1 0 0 0]")
    write("p", repr(a.p), "[1 -1 -2 0 0 0 0]")
    # Placeholder; recalculated from p and T (calculateDensity yes)
    write("rho", "1", "[1 -3 0 0 0 0 0]")
    write("U", "(0 0 0)", "[0 1 -1 0 0 0 0]", "volVectorField", "slip")
    import re
    cd = open(os.path.join("system", "controlDict")).read()
    cd = re.sub(r"\ndeltaT\s+[^;]+;", f"\ndeltaT          {a.deltaT!r};", cd)
    cd = re.sub(r"\nendTime\s+[^;]+;", f"\nendTime         {a.endTime!r};", cd)
    cd = re.sub(r"\nwriteInterval\s+[^;]+;", f"\nwriteInterval   {a.endTime!r};", cd)
    open(os.path.join("system", "controlDict"), "w").write(cd)

    for k, y in mixture(a.Y).items():
        write(k, repr(y), "[0 0 0 0 0 0 0]")
    write("Ydefault", "0", "[0 0 0 0 0 0 0]")


if __name__ == "__main__":
    main()

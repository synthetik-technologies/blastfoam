#!/usr/bin/env python3
"""
setup_case.py: initial fields of the 0D equilibrium chemistry case
(uniform T, p, zero velocity, mass fractions of the chosen mixture).

    python3 setup_case.py --T 4000 --p 1e7 --mixture productsAir
"""

import argparse
import os

SPECIES = ["O2", "N2", "CO", "H2O", "CO2", "NH3", "CH4", "C2H6", "Al", "H2",
           "O", "N", "NO", "H", "OH", "AlO", "C(gr)", "Al2O3(s)"]
AIR = {"N2": 0.7669, "O2": 0.2331}
# Comp-B/Al detonation products (pressurePhaseChange target fractions)
PRODUCTS = {"H2": 0.00433265, "N2": 0.28739527, "CO": 0.28672579,
            "H2O": 0.11160467, "CO2": 0.15365548, "NH3": 0.00144386,
            "CH4": 0.02801544, "C(gr)": 0.04906294, "Al2O3(s)": 0.07776391}


def mixture(name):
    if name == "air":
        Y = dict(AIR)
    elif name == "products":
        Y = dict(PRODUCTS)
    elif name == "productsAir":
        s = sum(PRODUCTS.values())
        Y = {k: 0.5*v/s for k, v in PRODUCTS.items()}
        for k, v in AIR.items():
            Y[k] = Y.get(k, 0) + 0.5*v
    else:
        raise SystemExit(f"Unknown mixture {name}")
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
    ap.add_argument("--mixture", default="productsAir")
    a = ap.parse_args()

    os.makedirs("0", exist_ok=True)
    write("T", repr(a.T), "[0 0 0 1 0 0 0]")
    write("p", repr(a.p), "[1 -1 -2 0 0 0 0]")
    # Placeholder; recalculated from p and T (calculateDensity yes)
    write("rho", "1", "[1 -3 0 0 0 0 0]")
    write("U", "(0 0 0)", "[0 1 -1 0 0 0 0]", "volVectorField", "slip")
    for k, y in mixture(a.mixture).items():
        write(k, repr(y), "[0 0 0 0 0 0 0]")
    write("Ydefault", "0", "[0 0 0 0 0 0 0]")


if __name__ == "__main__":
    main()

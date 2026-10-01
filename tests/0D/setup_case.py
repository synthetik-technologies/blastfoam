#!/usr/bin/env python3
"""
setup_case.py: Generates initial condition fields for tests/0D/phaseTransfer_0D.
"""

import os

CASE_DIR = "/opt/OpenFOAM/blastfoam/tests/0D/phaseTransfer_0D"
ZERO_DIR = os.path.join(CASE_DIR, "0")
os.makedirs(ZERO_DIR, exist_ok=True)

HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  12
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    location    "0";
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      {dims};

internalField   uniform {val};

boundaryField
{{
    sides
    {{
        type            {b_sides};
    }}
    bottomAndTop
    {{
        type            empty;
    }}
    frontAndBack
    {{
        type            empty;
    }}
}}

// ************************************************************************* //
"""


def write_field(filename, obj, cls, dims, val, b_sides="symmetry"):
    filepath = os.path.join(ZERO_DIR, filename)
    with open(filepath, "w") as f:
        f.write(HEADER.format(
            cls=cls,
            obj=obj,
            dims=dims,
            val=val,
            b_sides=b_sides
        ))


# Write fields
write_field("p", "p", "volScalarField", "[1 -1 -2 0 0 0 0]", "1.0e7", "symmetry")
write_field("T", "T", "volScalarField", "[0 0 0 1 0 0 0]", "1318.69", "symmetry")
write_field("U", "U", "volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)", "symmetry")

write_field("alpha.compBAl", "alpha.compBAl", "volScalarField", "[0 0 0 0 0 0 0]", "1", "symmetry")
write_field("alpha.gas", "alpha.gas", "volScalarField", "[0 0 0 0 0 0 0]", "0", "symmetry")
write_field("lambda.compBAl", "lambda.compBAl", "volScalarField", "[0 0 0 0 0 0 0]", "1", "symmetry")

write_field("rho.compBAl", "rho.compBAl", "volScalarField", "[1 -3 0 0 0 0 0]", "38.137262", "symmetry")
write_field("rho.gas", "rho.gas", "volScalarField", "[1 -3 0 0 0 0 0]", "38.137262", "symmetry")

# Specific internal energy consistent with 100 bar transfer point (p=100 bar, e=430110.0 J/kg)
write_field("e", "e", "volScalarField", "[0 2 -2 0 0 0 0]", "430110.0", "symmetry")

# Gas species (EuDetCode freeze-out composition normalized to 1.0)
species_fractions = {
    "N2.gas": "0.28466454",
    "CO.gas": "0.23116680",
    "H2O.gas": "0.10059566",
    "CO2.gas": "0.21098004",
    "NH3.gas": "0.00548327",
    "CH4.gas": "0.04230543",
    "C2H6.gas": "0.00149695",
    "C(gr).gas": "0.04548207",
    "Al2O3(s).gas": "0.07782524",
    "O2.gas": "0.0"
}

for sp, val in species_fractions.items():
    write_field(sp, sp, "volScalarField", "[0 0 0 0 0 0 0]", val, "symmetry")

print("Generated 0/ directory fields in tests/0D/phaseTransfer_0D")

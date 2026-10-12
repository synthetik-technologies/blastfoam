#!/usr/bin/env python3
"""
setup_case.py: Generates initial condition fields for tests/0D/phaseTransfer_0D
(2-phase: compBAl products + gas) or, with --three-phase, for
tests/0D/phaseTransfer_0D_3phase (compBAl products + inert Al + gas).
"""

import os
import sys

THREE_PHASE = "--three-phase" in sys.argv[1:]

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
CASE_DIR = os.path.join(
    TESTS_DIR, "phaseTransfer_0D_3phase" if THREE_PHASE else "phaseTransfer_0D"
)
ZERO_DIR = os.path.join(CASE_DIR, "0")
os.makedirs(ZERO_DIR, exist_ok=True)

# 100 bar freeze-out state of the reactive products (eJWLC isentrope).
# T0 is the gas-EOS temperature at this state; the eJWLC omegaT is set so
# the products isentrope temperature gives the same value.
P0 = 1.0e7
T0 = 825.76
RHO_PROD = 38.137262

# Inert Al: 90% of the 30 wt% Al in 70/30 Comp-B/Al, i.e. 27% of total mass,
# in thermal and mechanical equilibrium with the products. Murnaghan and
# eConst parameters must match constant/physicalProperties.Al of the 3-phase case.
Y_AL = 0.27
AL_RHO0, AL_K0, AL_N, AL_PREF = 2700.0, 72.7e9, 4.3, 101325.0

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


if THREE_PHASE:
    # Murnaghan: p = pRef + K0/n*((rho/rho0)^n - 1)
    rho_al = AL_RHO0*(1.0 + AL_N*(P0 - AL_PREF)/AL_K0)**(1.0/AL_N)
    v_prod = (1.0 - Y_AL)/RHO_PROD
    v_al = Y_AL/rho_al
    alpha_al = v_al/(v_prod + v_al)
    alpha_comp = 1.0 - alpha_al
else:
    alpha_comp = 1.0

# Write fields (blastEulerFoam: per-phase alpha, rho, U and T; each phase
# energy is initialised from its own temperature)
phases = (["Al"] if THREE_PHASE else []) + ["compBAl", "gas"]
write_field("p", "p", "volScalarField", "[1 -1 -2 0 0 0 0]", "1.0e7", "symmetry")
write_field("T", "T", "volScalarField", "[0 0 0 1 0 0 0]", f"{T0}", "symmetry")
write_field("U", "U", "volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)", "symmetry")
for ph in phases:
    write_field(f"T.{ph}", f"T.{ph}", "volScalarField", "[0 0 0 1 0 0 0]", f"{T0}", "symmetry")
    write_field(f"U.{ph}", f"U.{ph}", "volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)", "symmetry")

write_field("alpha.compBAl", "alpha.compBAl", "volScalarField", "[0 0 0 0 0 0 0]", f"{alpha_comp:.10g}", "symmetry")
write_field("alpha.gas", "alpha.gas", "volScalarField", "[0 0 0 0 0 0 0]", "0", "symmetry")

write_field("rho.compBAl", "rho.compBAl", "volScalarField", "[1 -3 0 0 0 0 0]", f"{RHO_PROD}", "symmetry")
write_field("rho.gas", "rho.gas", "volScalarField", "[1 -3 0 0 0 0 0]", f"{RHO_PROD}", "symmetry")

if THREE_PHASE:
    write_field("alpha.Al", "alpha.Al", "volScalarField", "[0 0 0 0 0 0 0]", f"{alpha_al:.10g}", "symmetry")
    write_field("rho.Al", "rho.Al", "volScalarField", "[1 -3 0 0 0 0 0]", f"{rho_al:.8f}", "symmetry")
    write_field("Theta.Al", "Theta.Al", "volScalarField", "[0 2 -2 0 0 0 0]", "0", "symmetry")

write_field("Ydefault.gas", "Ydefault.gas", "volScalarField", "[0 0 0 0 0 0 0]", "0", "symmetry")

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

print(f"Generated 0/ directory fields in {os.path.relpath(CASE_DIR, os.path.dirname(TESTS_DIR))}")

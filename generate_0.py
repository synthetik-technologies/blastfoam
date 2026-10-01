import os

def write_foam_file(path, name, cls, dims, internal, bc_type="zeroGradient"):
    content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  12
     \\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    location    "0";
    object      {name};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      {dims};

internalField   {internal};

boundaryField
{{
    inlet
    {{
        type            {bc_type};
    }}
    outlet
    {{
        type            {bc_type};
    }}
    bottom
    {{
        type            empty;
    }}
    top
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
    with open(path, 'w') as f:
        f.write(content)

os.makedirs('tutorials/blastFoam/multiphase_reacting_bubble_1d/0', exist_ok=True)
base = 'tutorials/blastFoam/multiphase_reacting_bubble_1d/0/'

write_foam_file(base + 'alpha.air', 'alpha.air', 'volScalarField', '[0 0 0 0 0 0 0]', 'uniform 1')
write_foam_file(base + 'alpha.bubble', 'alpha.bubble', 'volScalarField', '[0 0 0 0 0 0 0]', 'uniform 0')
write_foam_file(base + 'CH4.bubble', 'CH4.bubble', 'volScalarField', '[0 0 0 0 0 0 0]', 'uniform 0')
write_foam_file(base + 'CO2.bubble', 'CO2.bubble', 'volScalarField', '[0 0 0 0 0 0 0]', 'uniform 0')
write_foam_file(base + 'H2O.bubble', 'H2O.bubble', 'volScalarField', '[0 0 0 0 0 0 0]', 'uniform 0')
write_foam_file(base + 'N2.bubble', 'N2.bubble', 'volScalarField', '[0 0 0 0 0 0 0]', 'uniform 1')
write_foam_file(base + 'O2.bubble', 'O2.bubble', 'volScalarField', '[0 0 0 0 0 0 0]', 'uniform 0')
write_foam_file(base + 'p', 'p', 'volScalarField', '[1 -1 -2 0 0 0 0]', 'uniform 1e5')
write_foam_file(base + 'rho.air', 'rho.air', 'volScalarField', '[1 -3 0 0 0 0 0]', 'uniform 0.1742')
write_foam_file(base + 'rho.bubble', 'rho.bubble', 'volScalarField', '[1 -3 0 0 0 0 0]', 'uniform 0.1732')
write_foam_file(base + 'T', 'T', 'volScalarField', '[0 0 0 1 0 0 0]', 'uniform 2000')
write_foam_file(base + 'U', 'U', 'volVectorField', '[0 1 -1 0 0 0 0]', 'uniform (100 0 0)')

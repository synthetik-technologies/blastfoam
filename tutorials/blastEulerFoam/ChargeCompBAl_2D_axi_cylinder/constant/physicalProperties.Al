/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  12
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      physicalProperties.Al;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

type basic;

thermoType
{
    equationOfState Murnaghan;
    thermo          eMelting;
    transport       constIso;
}

specie
{
    molWeight       26.98;
}

equationOfState
{
    rho0            2700;
    pRef            101325;
    K0              72.7e9;
    n               4.3;
    Gamma           2.0;
}

transport
{
    kappa           237;
}

thermodynamics
{
    // NASA (CODATA) Al(cr)/Al(L), tests/0D/alMelting_0D
    type            eMelting;
    Cv              1059;
    CvLiquid        1177;
    Tm              933.47;
    Lm              3.966e5;
    dTm             2;
    Hf              0;
}

residualRho     2700;
residualAlpha   1e-06;

// ************************************************************************* //

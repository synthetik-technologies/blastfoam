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
    thermo          eConst;
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
    type            eConst;
    Cv              900;
    Hf              0;
}

residualRho     2700;
residualAlpha   1e-09;

// ************************************************************************* //

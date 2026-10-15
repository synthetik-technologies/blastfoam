/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  12
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      rho.Al;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [1 -3 0 0 0 0 0];

internalField   uniform 2700.181911;

boundaryField
{
    sides
    {
        type            symmetry;
    }
    bottomAndTop
    {
        type            empty;
    }
    frontAndBack
    {
        type            empty;
    }
}

// ************************************************************************* //

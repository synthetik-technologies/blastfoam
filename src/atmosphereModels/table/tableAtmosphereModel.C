/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technologies
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "tableAtmosphereModel.H"
#include "thermodynamicConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace atmosphereModels
{
    defineTypeNameAndDebug(table, 0);
    addToRunTimeSelectionTable(atmosphereModel, table, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atmosphereModels::table::table
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    atmosphereModel(mesh, dict),
    pTable_(dict_.subDict("pTable"), "h", "p"),
    TTable_(dict_.subDict("TTable"), "h", "T"),
    correct_(dict_.lookupOrDefault("correct", false))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::atmosphereModels::table::~table()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::atmosphereModels::table::createAtmosphere
(
    fluidBlastThermo& thermo
) const
{
    volScalarField& p = thermo.p();
    volScalarField& T = thermo.T();

    forAll(p, celli)
    {
        p[celli] = pTable_.lookup(h_[celli]);
        T[celli] = TTable_.lookup(h_[celli]);
    }

    p.correctBoundaryConditions();
    T.correctBoundaryConditions();

    volScalarField::Boundary& bp = p.boundaryFieldRef();
    volScalarField::Boundary& bT = T.boundaryFieldRef();
    forAll(bp, patchi)
    {
        if (bp[patchi].fixesValue())
        {
            forAll(bp[patchi], facei)
            {
                bp[patchi][facei] =
                    pTable_.lookup(h_.boundaryField()[patchi][facei]);
            }
        }
        if (bT[patchi].fixesValue())
        {
            forAll(bT[patchi], facei)
            {
                bT[patchi][facei] =
                    TTable_.lookup(h_.boundaryField()[patchi][facei]);
            }
        }
    }

    // Correct density
    thermo.updateRho(p);

    // Equalibriate the pressure field
    if (correct_)
    {
        hydrostaticInitialisation
        (
            thermo,
            dimensionedScalar(dimPressure, pTable_.lookup(gMin(h_)))
        );
    }
}

// ************************************************************************* //

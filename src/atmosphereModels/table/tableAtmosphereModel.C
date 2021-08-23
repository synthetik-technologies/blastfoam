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
    pTable_
    (
        dict_.subDict("pTable").lookup<fileName>("file"),
        dict_.subDict("pTable").lookupOrDefault<word>("pMod", "none"),
        dict_.subDict("pTable").lookupOrDefault<word>("hMod", "none")
    ),
    TTable_
    (
        dict_.subDict("TTable").lookup<fileName>("file"),
        dict_.subDict("TTable").lookupOrDefault<word>("TMod", "none"),
        dict_.subDict("TTable").lookupOrDefault<word>("hMod", "none")
    ),
    W_(dict_.lookupOrDefault("W", 28.97))
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
    forAll(thermo.p(), celli)
    {
        thermo.p()[celli] = pTable_.lookup(h_[celli]);
        thermo.T()[celli] = TTable_.lookup(h_[celli]);
    }

    thermo.calce(thermo.p());
    thermo.updateRho();
    hydrostaticInitialisation
    (
        thermo,
        dimensionedScalar(dimPressure, pTable_.lookup(gMin(h_)))
    );
}

// ************************************************************************* //

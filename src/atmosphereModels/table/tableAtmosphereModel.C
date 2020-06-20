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
        dict_.subDict("pTable").lookup("file"),
        dict_.subDict("pTable").lookupOrDefault<word>("pMod", "none"),
        dict_.subDict("pTable").lookupOrDefault<word>("hMod", "none")
    ),
    TTable_
    (
        dict_.subDict("TTable").lookup("file"),
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
    volScalarField& p,
    volScalarField& rho
) const
{
    scalar R = Foam::constant::thermodynamic::RR/W_;
    forAll(rho, celli)
    {
        p[celli] = pTable_.lookup(h_[celli]);
        scalar T = TTable_.lookup(h_[celli]);
        rho[celli] = p[celli]/T/R;
    }

    forAll(rho.boundaryField(), patchi)
    {
        forAll(rho.boundaryField()[patchi], facei)
        {
            p.boundaryFieldRef()[patchi][facei] =
                pTable_.lookup(h_.boundaryField()[patchi][facei]);
            scalar T = TTable_.lookup(h_.boundaryField()[patchi][facei]);
            rho.boundaryFieldRef()[patchi][facei] =
                p.boundaryField()[patchi][facei]/T/R;
        }
    }
}

// ************************************************************************* //

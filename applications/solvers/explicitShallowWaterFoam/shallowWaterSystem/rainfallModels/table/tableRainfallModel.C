/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "tableRainfallModel.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace rainfallModels
{
    defineTypeNameAndDebug(table, 0);
    addToRunTimeSelectionTable(rainfallModel, table, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rainfallModels::table::table
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    rainfallModel(mesh, dict),
    table_(dict, "x", "y", "t", "R")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rainfallModels::table::~table()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::rainfallModels::table::R0() const
{
    tmp<volScalarField> tR
    (
        volScalarField::New
        (
            typeName + ":R0",
            mesh_,
            dimensionedScalar(dimLength/dimTime, 0.0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    const scalar t = mesh_.time().value();
    if (t < start_ || t > end_)
    {
        return tR;
    }

    volScalarField& R = tR.ref();
    const volVectorField& C = mesh_.C();
    forAll(R, celli)
    {
        R[celli] = table_.lookup(C[celli].x(), C[celli].y(), t);
    }
    R.correctBoundaryConditions();
    return tR;
}

// ************************************************************************* //

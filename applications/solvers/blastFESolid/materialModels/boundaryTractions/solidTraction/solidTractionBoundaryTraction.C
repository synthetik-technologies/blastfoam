/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
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

#include "solidTractionBoundaryTraction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace boundaryTractions
{
    defineTypeNameAndDebug(solidTraction, 0);
    addToRunTimeSelectionTable
    (
        boundaryTraction,
        solidTraction,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryTractions::solidTraction::solidTraction
(
    const word& name,
    const dictionary& dict,
    const feMesh1& femesh,
    const pointVectorField* DPtr
)
:
    boundaryTraction(name, dict, femesh, DPtr),
    pressure_
    (
        Function1<scalar>::New
        (
            "pressure",
            femesh.time().userUnits(),
            dimPressure,
            dict
        )
    ),
    traction_
    (
        Function1<vector>::New
        (
            "traction",
            femesh.time().userUnits(),
            dimPressure,
            dict
        )
    )
{}


Foam::boundaryTractions::solidTraction::solidTraction
(
    const solidTraction& stbt
)
:
    boundaryTraction(stbt),
    pressure_(stbt.pressure_, false),
    traction_(stbt.traction_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boundaryTractions::solidTraction::~solidTraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::boundaryTractions::solidTraction::pressure
(
    const labelList& nodes,
    const scalarList& shape,
    const vector& n
) const
{
   return pressure_->value(mesh_.time().value());
}


Foam::vector Foam::boundaryTractions::solidTraction::traction
(
    const labelList& nodes,
    const scalarList& shape,
    const vector& n
) const
{
   return traction_->value(mesh_.time().value());
}


// ************************************************************************* //

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

#include "fixedVelocityDisplacementConstraint.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace displacementConstraints
{
    defineTypeNameAndDebug(fixedVelocity, 0);
    addToRunTimeSelectionTable
    (
        displacementConstraint,
        fixedVelocity,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementConstraints::fixedVelocity::fixedVelocity
(
    const word& name,
    pointVectorField& D,
    pointVectorField& U,
    const dictionary& dict
)
:
    displacementConstraint(name, D, U, dict),
    velocity_
    (
        Function1<vector>::New
        (
            "velocity",
            D.time().userUnits(),
            dimVelocity,
            dict
        )
    )
{}


Foam::displacementConstraints::fixedVelocity::fixedVelocity
(
    const fixedVelocity& fvdc
)
:
    displacementConstraint(fvdc),
    velocity_(fvdc.velocity_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementConstraints::fixedVelocity::~fixedVelocity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::displacementConstraints::fixedVelocity::constrain()
{
    const scalar t = D_.time().value();
    const scalar dt = D_.time().deltaTValue();
    const pointVectorField& D0 = D_.oldTime();

    const vector u = velocity_->value(t);
    forAll(nodes_, ni)
    {
        const label nodei = nodes_[ni];
        forAll(dims_, ci)
        {
            const label cmpti = dims_[ci];
            D_[nodei][cmpti] = D0[nodei][cmpti] + dt*u[cmpti];
            U_[nodei][cmpti] = u[cmpti];
        }
    }
}

// ************************************************************************* //

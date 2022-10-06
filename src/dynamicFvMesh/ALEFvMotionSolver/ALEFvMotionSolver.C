/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020   21-05-2020 Synthetik Applied Technologies: |  Set old cell volumes
                                after refinement is done.
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

#include "ALEFvMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ALEFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        ALEFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ALEFvMotionSolver::ALEFvMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    motionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ALEFvMotionSolver::~ALEFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::ALEFvMotionSolver::curPoints() const
{
    const volVectorField& U = fvMesh_.lookupObject<volVectorField>("U");
    pointVectorField pU
    (
        pointVectorField::New
        (
            "pU",
            pointMesh::New(fvMesh_),
            dimensionedVector(dimVelocity, Zero),
            wordList(U.boundaryField().size(), "slip")
        )
    );
    volPointInterpolation::New(fvMesh_).interpolate(U, pU);
    pU.correctBoundaryConditions();
    tmp<pointField> tcurPoints
    (
        fvMesh_.points()
      + fvMesh_.time().deltaTValue()*pU.primitiveField()
    );

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


// ************************************************************************* //

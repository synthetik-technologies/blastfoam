/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
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

#include "noAfterburn.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace afterburnModels
{
    defineTypeNameAndDebug(noAfterburn, 0);
    addToRunTimeSelectionTable(afterburnModel, noAfterburn, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::afterburnModels::noAfterburn::noAfterburn
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    afterburnModel(mesh, dict, phaseName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::afterburnModels::noAfterburn::~noAfterburn()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField>
Foam::afterburnModels::noAfterburn::ESource() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "noAfterburn:ESource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("0", sqr(dimVelocity)/dimTime, 0.0)
        )
    );
}

// ************************************************************************* //

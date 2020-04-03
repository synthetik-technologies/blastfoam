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

#include "linearAfterburn.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace afterburnModels
{
    defineTypeNameAndDebug(linearAfterburn, 0);
    addToRunTimeSelectionTable(afterburnModel, linearAfterburn, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::afterburnModels::linearAfterburn::linearAfterburn
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    afterburnModel(mesh, dict, phaseName),
    Q0_("Q0", sqr(dimVelocity), dict),
    tStart_(dimensionedScalar::lookupOrDefault("tStart", dict, dimTime, 0.0)),
    tEnd_("tEnd", dimTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::afterburnModels::linearAfterburn::~linearAfterburn()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::afterburnModels::linearAfterburn::ESource() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "linearAfterburn:Esource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            Q0_/(tEnd_ - tStart_)
           *pos(time_ - tStart_)*pos(tEnd_ - tStart_)
        )
    );
}

// ************************************************************************* //

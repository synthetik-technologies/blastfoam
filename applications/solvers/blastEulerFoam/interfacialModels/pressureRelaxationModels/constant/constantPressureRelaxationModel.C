/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "constantPressureRelaxationModel.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pressureRelaxationModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable(pressureRelaxationModel, constant, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureRelaxationModels::constant::constant
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    pressureRelaxationModel(dict, pair, registerObject),
    mu_("mu", dimensionSet(-1, 1, 1, 0, 0, 0, 0), dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pressureRelaxationModels::constant::~constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::pressureRelaxationModels::constant::K
(
    const label nodei,
    const label nodej
) const
{
    const fvMesh& mesh = pair_.phase1().mesh();
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "constant:K",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            mu_
        )
    );
}


Foam::scalar Foam::pressureRelaxationModels::constant::Ki
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    return mu_.value();
}

// ************************************************************************* //

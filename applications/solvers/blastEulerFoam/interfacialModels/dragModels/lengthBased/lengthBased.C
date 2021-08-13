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

#include "lengthBased.H"
#include "phasePair.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(lengthBased, 0);
    addToRunTimeSelectionTable(dragModel, lengthBased, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::lengthBased::lengthBased
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    C_("C", dimArea/dimTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::lengthBased::~lengthBased()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::lengthBased::CdRe
(
    const label nodei,
    const label nodej
) const
{
    FatalErrorInFunction
        << "Not implemented."
        << "Drag coefficient not defined for the lengthBased model."
        << exit(FatalError);

    return pair_.phase1();
}


Foam::tmp<Foam::volScalarField> Foam::dragModels::lengthBased::KI
(
    const label nodei,
    const label nodej
) const
{
    return
        C_
       *pair_.dispersed().rho()
       *4.0
       /sqr(pair_.dispersed().d(nodei));
}


Foam::scalar Foam::dragModels::lengthBased::CdRei
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    FatalErrorInFunction
        << "Not implemented."
        << "Drag coefficient not defined for the lengthBased model."
        << exit(FatalError);

    return pair_.phase1()[celli];
}


Foam::scalar Foam::dragModels::lengthBased::KIi
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    return
        C_.value()
       *pair_.dispersed().rho()[celli]
       *4.0
       /sqr(pair_.dispersed().di(celli, nodei));
}

// ************************************************************************* //

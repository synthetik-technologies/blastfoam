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

#include "Stokes.H"
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
    defineTypeNameAndDebug(Stokes, 0);
    addToRunTimeSelectionTable(dragModel, Stokes, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Stokes::Stokes
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    dragTime_("dragTime", dimTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::Stokes::~Stokes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::Stokes::CdRe
(
    const label nodei,
    const label nodej
) const
{
    FatalErrorInFunction
        << "Not implemented."
        << "Drag coefficient not defined for the Stokes model."
        << exit(FatalError);

    return pair_.phase1();
}


Foam::tmp<Foam::volScalarField> Foam::dragModels::Stokes::Ki
(
    const label nodei,
    const label nodej
) const
{
    return pair_.dispersed().rho()/dragTime_;
}


Foam::scalar Foam::dragModels::Stokes::CdRe
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    FatalErrorInFunction
        << "Not implemented."
        << "Drag coefficient not defined for the Stokes model."
        << exit(FatalError);

    return pair_.phase1()[celli];
}


Foam::scalar Foam::dragModels::Stokes::Ki
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    return pair_.dispersed().rho()[celli]/dragTime_.value();
}


// ************************************************************************* //

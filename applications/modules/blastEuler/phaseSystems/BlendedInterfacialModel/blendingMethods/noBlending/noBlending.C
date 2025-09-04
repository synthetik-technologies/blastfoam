/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "noBlending.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(noBlending, 0);

    addToRunTimeSelectionTable
    (
        blendingMethod,
        noBlending,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::noBlending::noBlending
(
    const dictionary& dict,
    const wordList& phaseNames
)
:
    blendingMethod(dict),
    continuousPhase_(dict.lookupOrDefault("continuousPhase", word::null)),
    dispersedPhase_(dict.lookupOrDefault("dispersedPhase", word::null))
{
    if (continuousPhase_.empty() && dispersedPhase_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Either continuousPhase or dispersedPhase must be provided"
            << endl
            << abort(FatalIOError);
    }
    else if (!continuousPhase_.empty() && !dispersedPhase_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Either continuousPhase or dispersedPhase must be provided, "
            << "not both" << endl
            << abort(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::noBlending::~noBlending()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::noBlending::f1
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const fvMesh& mesh(phase1.mesh());

    return volScalarField::New
    (
        "f",
        mesh,
        continuousPhase_.empty()
      ? dimensionedScalar(dimless, phase1.name() == dispersedPhase_)
      : dimensionedScalar(dimless, phase2.name() == continuousPhase_)
    );
}


Foam::scalar Foam::blendingMethods::noBlending::f1i
(
    const label celli,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return scalar
    (
        continuousPhase_.empty()
      ? phase1.name() == dispersedPhase_
      : phase2.name() == continuousPhase_
    );
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::noBlending::f2
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const fvMesh& mesh(phase1.mesh());

    return volScalarField::New
    (
        "f",
        mesh,
        continuousPhase_.empty()
      ? dimensionedScalar(dimless, phase2.name() == dispersedPhase_)
      : dimensionedScalar(dimless, phase1.name() == continuousPhase_)
    );
}

Foam::scalar Foam::blendingMethods::noBlending::f2i
(
    const label celli,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return scalar
    (
        continuousPhase_.empty()
      ? phase2.name() == dispersedPhase_
      : phase1.name() == continuousPhase_
    );
}


// ************************************************************************* //

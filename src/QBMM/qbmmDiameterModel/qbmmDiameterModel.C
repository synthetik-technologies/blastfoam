/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technologies
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

#include "qbmmDiameterModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(qbmmDiameterModel, 0);
    addToRunTimeSelectionTable(diameterModel, qbmmDiameterModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::qbmmDiameterModel::qbmmDiameterModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    diameterModel(mesh, dict, phaseName),
    pbeDict_
    (
        IOobject
        (
            IOobject::groupName("populationBalanceProperties", phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ
        )
    ),
    pbe_
    (
        phaseName,
        pbeDict_.subDict("univariateCoeffs"),
        mesh.foundObject<surfaceScalarField>
        (
            IOobject::groupName("phi", phaseName)
        )
        ? mesh.lookupObject<surfaceScalarField>
        (
            IOobject::groupName("phi", phaseName)
        )
        : mesh.lookupObject<surfaceScalarField>("phi")
    ),
    numberDensity_
    (
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("moment.0", phaseName)
        ).dimensions() == inv(dimVolume)
    ),
    vf_
    (
        numberDensity_
      ? mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("moment.3", phaseName)
        )
      : mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("moment.0", phaseName)
        )
    ),
    vfD_
    (
        numberDensity_
      ? mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("moment.4", phaseName)
        )
      : mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("moment.1", phaseName)
        )
    )
{
    if (vfD_.dimensions() != dimLength)
    {
        FatalErrorInFunction
            << "The primary abscissa must describe the diameter of the" << nl
            << "dispersed phase." << endl
            << abort(FatalError);
    }
    this->d_ = vfD_/max(vf_, 1e-10);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::qbmmDiameterModel::~qbmmDiameterModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::qbmmDiameterModel::setODEFields
(
    const label nSteps,
    const labelList& oldIs,
    const label& nOld,
    const labelList& deltaIs,
    const label nDelta
)
{
    boolList storeOld(nSteps, false);
    forAll(oldIs, i)
    {
        if (oldIs[i] >= 0)
        {
            storeOld[i] = true;
        }
    }
    boolList storeDeltas(nSteps, false);
    forAll(deltaIs, i)
    {
        if (deltaIs[i] >= 0)
        {
            storeDeltas[i] = true;
        }
    }

    pbe_.setODEFields(nSteps, storeOld, storeDeltas);
}


void Foam::diameterModels::qbmmDiameterModel::clearODEFields()
{
    pbe_.clearODEFields();
}


void Foam::diameterModels::qbmmDiameterModel::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    pbe_.solve(stepi, ai, bi);
    this->d_ = vfD_/max(vf_, 1e-10);
}

// ************************************************************************* //

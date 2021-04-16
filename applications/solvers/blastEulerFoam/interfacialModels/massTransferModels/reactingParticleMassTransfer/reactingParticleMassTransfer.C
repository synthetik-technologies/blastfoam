/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "reactingParticleMassTransfer.H"
#include "phaseSystem.H"
#include "phasePair.H"
#include "interfacialPressureModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferModels
{
    defineTypeNameAndDebug(reactingParticleMassTransfer, 0);
    addToRunTimeSelectionTable(massTransferModel, reactingParticleMassTransfer, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massTransferModels::reactingParticleMassTransfer::
reactingParticleMassTransfer
(
    const dictionary& dict,
    const phasePair& pair
)
:
    massTransferModel(dict, pair),
    reactingSpecies_(dict.lookupOrDefault("reactingSpecies", wordList())),
    productSpecies_(dict.lookupOrDefault("productSpecies", wordList()))
{
    scalarList rY;
    scalarList pY;
    if (reactingSpecies_.size())
    {
        rY = scalarList(dict.lookup("reactingYi"));
    }
    if (productSpecies_.size())
    {
        pY = scalarList(dict.lookup("productYi"));
    }

    forAll(reactingSpecies_, i)
    {
        reactingYi_.insert(reactingSpecies_[i], rY[i]);
    }
    forAll(productSpecies_, i)
    {
        productYi_.insert(productSpecies_[i], pY[i]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::massTransferModels::reactingParticleMassTransfer::~reactingParticleMassTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::massTransferModels::reactingParticleMassTransfer::K() const
{
    const diameterModel& dModel = pair_.dispersed().dModel();
    volScalarField V(dModel.V());
    V.max(small);
    tmp<volScalarField> n(pair_.dispersed()/V);

    tmp<volScalarField> mDot
    (
        volScalarField::New
        (
            "reactingParticle:mDot",
            n
           *dModel.dMdt()
        )
    );
    return mDot;
}


Foam::tmp<Foam::volScalarField>
Foam::massTransferModels::reactingParticleMassTransfer::dispersedYi
(
    const word& name
) const
{
    scalar value = 0;
    if (reactingSpecies_.found(name))
    {
        value = reactingYi_[name];
    }
    return volScalarField::New
    (
        IOobject::groupName("Yi", name),
        pair_.phase1().mesh(),
        dimensionedScalar(dimless, value)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::massTransferModels::reactingParticleMassTransfer::continuousYi
(
    const word& name
) const
{
    scalar value = 0;
    if (productSpecies_.found(name))
    {
        value = productYi_[name];
    }
    return volScalarField::New
    (
        IOobject::groupName("Yi", name),
        pair_.phase1().mesh(),
        dimensionedScalar(dimless, value)
    );
}

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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
    massTransferModel(dict, pair)
{
    if (dict.found("reactingSpecies"))
    {
        List<Tuple2<word, scalar>> table(dict.lookup("reactingSpecies"));
        forAll(table, i)
        {
            reactingYi_.insert(table[i].first(), table[i].second());
            reactingSpecies_.append(table[i].first());
        }
    }
    if (dict.found("productSpecies"))
    {
        List<Tuple2<word, scalar>> table(dict.lookup("productSpecies"));
        forAll(table, i)
        {
            productYi_.insert(table[i].first(), table[i].second());
            productSpecies_.append(table[i].first());
        }
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

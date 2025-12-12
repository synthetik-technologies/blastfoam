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
    massTransferModel(dict, pair),
    reactingPhase_
    (
        pair.ordered()
      ? pair.dispersed().name()
      : dict.lookup<word>("reactingPhase")
    ),
    reactingIs1_(pair.phase1().name() == reactingPhase_)
{
    const phaseModel& phase1 = pair.phase1();
    if (dict.found(phase1.name() + "Species"))
    {
        List<Tuple2<word, scalar>> table(dict.lookup(phase1.name() + "Species"));
        forAll(table, i)
        {
            phase1Yi_.insert(table[i].first(), table[i].second());
            phase1Species_.append(table[i].first());
        }
    }
    else if (dict.found("reactingSpecies"))
    {
        HashTable<scalar, word>& reactingYi =
            reactingIs1_ ? phase1Yi_ : phase2Yi_;
        hashedWordList& reactingSpecies =
            reactingIs1_ ? phase1Species_ : phase2Species_;

        List<Tuple2<word, scalar>> table(dict.lookup("reactingSpecies"));
        forAll(table, i)
        {
            reactingYi.insert(table[i].first(), table[i].second());
            reactingSpecies.append(table[i].first());
        }
    }

    const phaseModel& phase2 = pair.phase2();
    if (dict.found(phase2.name() + "Species"))
    {
        List<Tuple2<word, scalar>> table(dict.lookup(phase2.name() + "Species"));
        forAll(table, i)
        {
            phase2Yi_.insert(table[i].first(), table[i].second());
            phase2Species_.append(table[i].first());
        }
    }
    else if (dict.found("productSpecies"))
    {
        HashTable<scalar, word>& productYi =
            reactingIs1_ ? phase2Yi_ : phase1Yi_;
        hashedWordList& productSpecies =
            reactingIs1_ ? phase2Species_ : phase1Species_;

        List<Tuple2<word, scalar>> table(dict.lookup("productSpecies"));
        forAll(table, i)
        {
            productYi.insert(table[i].first(), table[i].second());
            productSpecies.append(table[i].first());
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
    return reactingPhase().dModel().dMdt();
}


Foam::tmp<Foam::volScalarField>
Foam::massTransferModels::reactingParticleMassTransfer::phase1Y
(
    const word& name
) const
{
    scalar value = 0;
    if (phase1Species_.found(name))
    {
        value = phase1Yi_[name];
    }
    return volScalarField::New
    (
        IOobject::groupName("Yi", name),
        pair_.phase1().mesh(),
        dimensionedScalar(dimless, value)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::massTransferModels::reactingParticleMassTransfer::phase2Y
(
    const word& name
) const
{
    scalar value = 0;
    if (phase2Species_.found(name))
    {
        value = phase2Yi_[name];
    }
    return volScalarField::New
    (
        IOobject::groupName("Yi", name),
        pair_.phase1().mesh(),
        dimensionedScalar(dimless, value)
    );
}

// ************************************************************************* //

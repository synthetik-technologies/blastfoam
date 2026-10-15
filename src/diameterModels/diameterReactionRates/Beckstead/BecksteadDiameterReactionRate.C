/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "BecksteadDiameterReactionRate.H"
#include "diameterModel.H"
#include "IOdictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterReactionRates
{
    defineTypeNameAndDebug(Beckstead, 0);
    addToRunTimeSelectionTable(diameterReactionRate, Beckstead, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::diameterReactionRates::Beckstead::readSpecies()
{
    const fvMesh& mesh = dModel_.d().mesh();
    const word dictName
    (
        IOobject::groupName("physicalProperties", oxidiserPhase_)
    );

    autoPtr<IOdictionary> readDict;
    const dictionary* propsPtr = nullptr;
    if (mesh.foundObject<IOdictionary>(dictName))
    {
        propsPtr = &mesh.lookupObject<IOdictionary>(dictName);
    }
    else
    {
        readDict.reset
        (
            new IOdictionary
            (
                IOobject
                (
                    dictName,
                    mesh.time().constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            )
        );
        propsPtr = &readDict();
    }
    const dictionary& props = *propsPtr;

    species_ = props.lookup<wordList>("species");
    W_.setSize(species_.size());
    forAll(species_, i)
    {
        W_[i] =
            props.subDict(species_[i]).subDict("specie").lookup<scalar>
            (
                "molWeight"
            );
    }

    // Oxidisers that are not species of the phase do not contribute
    bool found = false;
    forAllConstIter(HashTable<scalar>, oxidisers_, iter)
    {
        found = found || findIndex(species_, iter.key()) >= 0;
    }
    if (!found)
    {
        FatalErrorInFunction
            << "None of the oxidisers " << oxidisers_.toc()
            << " is a species of phase " << oxidiserPhase_
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterReactionRates::Beckstead::Beckstead
(
    const diameterModel& dModel,
    const dictionary& dict
)
:
    diameterReactionRate(dModel),
    a_(dict.lookupOrDefault<scalar>("a", 0.00735)),
    n_(dict.lookupOrDefault<scalar>("n", 1.8)),
    Tign_(dict.lookup<scalar>("Tign")),
    oxidiserPhase_(dict.lookup<word>("oxidiserPhase")),
    oxidisers_(),
    species_(),
    W_()
{
    if (dict.found("oxidisers"))
    {
        const dictionary& oxDict = dict.subDict("oxidisers");
        forAllConstIter(dictionary, oxDict, iter)
        {
            oxidisers_.insert(iter().keyword(), oxDict.lookup<scalar>(iter().keyword()));
        }
    }
    else
    {
        oxidisers_.insert("O2", 1.0);
        oxidisers_.insert("H2O", 0.6);
        oxidisers_.insert("CO2", 0.22);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterReactionRates::Beckstead::~Beckstead()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::diameterReactionRates::Beckstead::Xeff() const
{
    if (species_.empty())
    {
        const_cast<Beckstead&>(*this).readSpecies();
    }

    const fvMesh& mesh = dModel_.d().mesh();

    // Moles per unit mass of the oxidiser phase
    volScalarField moles
    (
        volScalarField::New
        (
            "moles",
            mesh,
            dimensionedScalar(dimless, small)
        )
    );
    volScalarField molesOx
    (
        volScalarField::New
        (
            "molesOx",
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );
    forAll(species_, i)
    {
        const volScalarField& Y =
            mesh.lookupObject<volScalarField>
            (
                IOobject::groupName(species_[i], oxidiserPhase_)
            );
        const volScalarField YbyW(max(Y, scalar(0))/W_[i]);
        moles += YbyW;
        if (oxidisers_.found(species_[i]))
        {
            molesOx += oxidisers_[species_[i]]*YbyW;
        }
    }

    return volScalarField::New("Xeff", molesOx/moles);
}


Foam::tmp<Foam::volScalarField>
Foam::diameterReactionRates::Beckstead::dDdt
(
    const volScalarField& p,
    const volScalarField& Tp
) const
{
    const volScalarField& d = dModel_.d();
    const fvMesh& mesh = d.mesh();
    const volScalarField& d0 =
        mesh.lookupObject<volScalarField>(IOobject::groupName("d0", d.group()));
    const volScalarField& Tgas =
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("T", oxidiserPhase_)
        );
    const volScalarField X(Xeff());

    tmp<volScalarField> tRate
    (
        volScalarField::New
        (
            IOobject::groupName("BecksteadRate", d.group()),
            mesh,
            dimensionedScalar(dimVelocity, 0)
        )
    );
    volScalarField& rate = tRate.ref();

    forAll(rate, celli)
    {
        const bool burning =
            Tp[celli] >= Tign_ || d[celli] < (1.0 - 1e-6)*d0[celli];

        if (burning && X[celli] > 0)
        {
            // dd/dt [um/ms] of the d^n law, converted to m/s; returned as
            // the rate of change of the radius
            const scalar dUm = max(d[celli], 1e-9)*1e6;
            const scalar dDdtUm =
                X[celli]*pow(max(p[celli], small)/101325.0, 0.1)
               *pow(max(Tgas[celli], small), 0.2)
               /(a_*n_*pow(dUm, n_ - 1.0));

            rate[celli] = 0.5*1e-3*dDdtUm;
        }
    }

    return tRate;
}


// ************************************************************************* //

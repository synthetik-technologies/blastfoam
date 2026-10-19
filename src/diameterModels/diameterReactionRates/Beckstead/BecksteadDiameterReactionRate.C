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

    if (!products_)
    {
        return;
    }

    // Mass yields of the reactions per unit mass of particle material. The
    // molecular weight of the particle material follows from the mass
    // balance of the reaction, so that the yields sum to one exactly.
    liquidOxide_ = findIndex(species_, oxideLiquid_) >= 0;
    forAllConstIter(HashTable<scalar>, oxidisers_, iter)
    {
        const word& ox = iter.key();
        if (findIndex(species_, ox) < 0)
        {
            continue;
        }
        if (!reactions_.found(ox))
        {
            FatalErrorInFunction
                << "No reaction for oxidiser " << ox << exit(FatalError);
        }
        const List<Tuple2<word, scalar>>& r = reactions_[ox];
        scalar mP = 0;
        forAll(r, i)
        {
            const label j = findIndex(species_, r[i].first());
            if (j < 0)
            {
                FatalErrorInFunction
                    << "Species " << r[i].first() << " of the " << ox
                    << " reaction is not a species of phase "
                    << oxidiserPhase_ << exit(FatalError);
            }
            mP += r[i].second()*W_[j];
        }
        const scalar Wp = mP/nuP_[ox];
        if (mag(Wp/26.981539 - 1) > 1e-3)
        {
            FatalErrorInFunction
                << "The " << ox << " reaction is not balanced for Al: "
                << "particle molecular weight " << Wp << exit(FatalError);
        }
        HashTable<scalar> y;
        forAll(r, i)
        {
            const label j = findIndex(species_, r[i].first());
            y.insert(r[i].first(), r[i].second()*W_[j]/mP);
            if (findIndex(productSpecies_, r[i].first()) < 0)
            {
                productSpecies_.append(r[i].first());
            }
        }
        yields_.insert(ox, y);
    }
    if
    (
        liquidOxide_
     && findIndex(productSpecies_, oxideSolid_) >= 0
     && findIndex(productSpecies_, oxideLiquid_) < 0
    )
    {
        productSpecies_.append(oxideLiquid_);
    }
}


void Foam::diameterReactionRates::Beckstead::shares
(
    PtrList<volScalarField>& w,
    wordList& names
) const
{
    if (species_.empty())
    {
        const_cast<Beckstead&>(*this).readSpecies();
    }

    const fvMesh& mesh = dModel_.d().mesh();
    names.clear();
    w.clear();
    volScalarField sum
    (
        volScalarField::New("sum", mesh, dimensionedScalar(dimless, 0))
    );
    forAllConstIter(HashTable<scalar>, oxidisers_, iter)
    {
        const label i = findIndex(species_, iter.key());
        if (i < 0)
        {
            continue;
        }
        const volScalarField& Y =
            mesh.lookupObject<volScalarField>
            (
                IOobject::groupName(species_[i], oxidiserPhase_)
            );
        names.append(iter.key());
        w.append
        (
            volScalarField::New
            (
                IOobject::groupName("share", iter.key()),
                iter()*max(Y, scalar(0))/W_[i]
            ).ptr()
        );
        sum += w.last();
    }
    forAll(w, k)
    {
        forAll(sum, celli)
        {
            w[k][celli] =
                sum[celli] > 0 ? w[k][celli]/sum[celli] : 1.0/w.size();
        }
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
    W_(),
    products_(dict.lookupOrDefault<bool>("products", false)),
    oxidiserLimit_(dict.lookupOrDefault<scalar>("oxidiserLimit", 0.5)),
    Textinction_(dict.lookupOrDefault<scalar>("Textinction", 0)),
    externalIgnition_(dict.lookupOrDefault<bool>("ignitedFraction", false)),
    oxideSolid_(dict.lookupOrDefault<word>("oxideSolid", "Al2O3(s)")),
    oxideLiquid_(dict.lookupOrDefault<word>("oxideLiquid", "Al2O3(L)")),
    TmOxide_(dict.lookupOrDefault<scalar>("TmOxide", 2327.0)),
    liquidOxide_(false)
{
    // Default global reactions (moles; particle Al, gas species)
    nuP_.insert("O2", 4);
    reactions_.insert
    (
        "O2",
        {Tuple2<word, scalar>("O2", -3), Tuple2<word, scalar>(oxideSolid_, 2)}
    );
    nuP_.insert("H2O", 2);
    reactions_.insert
    (
        "H2O",
        {
            Tuple2<word, scalar>("H2O", -3),
            Tuple2<word, scalar>(oxideSolid_, 1),
            Tuple2<word, scalar>("H2", 3)
        }
    );
    nuP_.insert("CO2", 2);
    reactions_.insert
    (
        "CO2",
        {
            Tuple2<word, scalar>("CO2", -3),
            Tuple2<word, scalar>(oxideSolid_, 1),
            Tuple2<word, scalar>("CO", 3)
        }
    );

    if (dict.found("oxidisers"))
    {
        const dictionary& oxDict = dict.subDict("oxidisers");
        forAllConstIter(dictionary, oxDict, iter)
        {
            const word& ox = iter().keyword();
            if (iter().isDict())
            {
                const dictionary& d = iter().dict();
                oxidisers_.insert(ox, d.lookup<scalar>("weight"));
                if (d.found("reaction"))
                {
                    nuP_.set(ox, d.lookup<scalar>("Al"));
                    reactions_.set
                    (
                        ox,
                        d.lookup<List<Tuple2<word, scalar>>>("reaction")
                    );
                }
            }
            else
            {
                oxidisers_.insert(ox, oxDict.lookup<scalar>(ox));
            }
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
            (
                externalIgnition_
             || Tp[celli] >= Tign_
             || d[celli] < (1.0 - 1e-6)*d0[celli]
            )
         && Tgas[celli] >= Textinction_;

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

    if (!products_)
    {
        return tRate;
    }

    // Limit the burning so that at most oxidiserLimit of each oxidiser in
    // the cell is consumed in a time step. The particle mass burnt per unit
    // volume and time is 6 alpha_p rho_p r'/d for a radius rate r'.
    PtrList<volScalarField> w;
    wordList oxNames;
    shares(w, oxNames);
    const volScalarField& alphaP =
        mesh.lookupObject<volScalarField>(IOobject::groupName("alpha", d.group()));
    const volScalarField& rhoP =
        mesh.lookupObject<volScalarField>(IOobject::groupName("rho", d.group()));
    const volScalarField& alphaG =
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", oxidiserPhase_)
        );
    const volScalarField& rhoG =
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("rho", oxidiserPhase_)
        );
    const scalar deltaT = mesh.time().deltaTValue();
    forAll(oxNames, k)
    {
        const word& ox = oxNames[k];
        const scalar mOx = -yields_[ox][ox];
        const volScalarField& Yox =
            mesh.lookupObject<volScalarField>
            (
                IOobject::groupName(ox, oxidiserPhase_)
            );
        forAll(rate, celli)
        {
            const scalar mBurnt =
                6*alphaP[celli]*rhoP[celli]*rate[celli]
               /max(d[celli], 1e-9);
            const scalar mOxRate = mBurnt*w[k][celli]*mOx;
            const scalar mOxMax =
                oxidiserLimit_*max(alphaG[celli], scalar(0))*rhoG[celli]
               *max(Yox[celli], scalar(0))/deltaT;
            if (mOxRate > mOxMax)
            {
                rate[celli] *= mOxMax/mOxRate;
            }
        }
    }

    return tRate;
}


Foam::wordList
Foam::diameterReactionRates::Beckstead::productSpecies() const
{
    if (!products_)
    {
        return wordList();
    }
    if (species_.empty())
    {
        const_cast<Beckstead&>(*this).readSpecies();
    }
    return productSpecies_;
}


Foam::tmp<Foam::volScalarField>
Foam::diameterReactionRates::Beckstead::productYield(const word& name) const
{
    tmp<volScalarField> tY(diameterReactionRate::productYield(name));
    if (!products_)
    {
        return tY;
    }
    volScalarField& Y = tY.ref();

    PtrList<volScalarField> w;
    wordList oxNames;
    shares(w, oxNames);

    // The oxide is solid or liquid according to the gas temperature
    const bool swap = liquidOxide_;
    const word yName =
        swap && name == oxideLiquid_ ? oxideSolid_ : name;
    forAll(oxNames, k)
    {
        const HashTable<scalar>& y = yields_[oxNames[k]];
        if (y.found(yName))
        {
            Y += y[yName]*w[k];
        }
        if (swap && name == oxideLiquid_ && y.found(oxideLiquid_))
        {
            // Liquid oxide named explicitly in a reaction
            Y += y[oxideLiquid_]*w[k];
        }
    }
    if (swap && (name == oxideSolid_ || name == oxideLiquid_))
    {
        const volScalarField& Tgas =
            dModel_.d().mesh().lookupObject<volScalarField>
            (
                IOobject::groupName("T", oxidiserPhase_)
            );
        forAll(Y, celli)
        {
            const bool liquid = Tgas[celli] >= TmOxide_;
            if (name == oxideSolid_ && liquid)
            {
                Y[celli] = 0;
            }
            else if (name == oxideLiquid_ && !liquid)
            {
                // Only the explicitly liquid products remain
                scalar yl = 0;
                forAll(oxNames, k)
                {
                    const HashTable<scalar>& y = yields_[oxNames[k]];
                    if (y.found(oxideLiquid_))
                    {
                        yl += y[oxideLiquid_]*w[k][celli];
                    }
                }
                Y[celli] = yl;
            }
        }
    }
    return tY;
}


// ************************************************************************* //

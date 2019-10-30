/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
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

#include "multiphaseEquationOfState.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseEquationOfState::multiphaseEquationOfState
(
    volScalarField& rho,
    volScalarField& e,
    volScalarField& p,
    const dictionary& dict
)
:
    phases_(dict.lookup("phases")),
    rho_(rho),
    e_(e),
    p_(p),
    eos_(phases_.size()),
    alphas_(phases_.size()),
    rhos_(phases_.size()),
    alphaRhos_(phases_.size()),
    alphaPhis_(phases_.size()),
    alphaRhoPhis_(phases_.size())
{
    forAll(phases_, phasei)
    {
        eos_.set
        (
            phasei,
            equationOfState::New(e_, dict.subDict(phases_[phasei])).ptr()
        );
        alphas_.set(phasei, &eos_[phasei]);
        rhos_.set(phasei, &eos_[phasei].rho());
        alphaRhos_.set(phasei, &eos_[phasei].alphaRho());
        alphaPhis_.set(phasei, &eos_[phasei].alphaPhi());
        alphaRhoPhis_.set(phasei, &eos_[phasei].alphaRhoPhi());

        rho_ += alphas_[phasei]*rhos_[phasei];
    }

    if (max(e_).value() < small)
    {
        forAll(phases_, phasei)
        {
            e_ += alphas_[phasei]*eos_[phasei].e(p_);
        }
        e_.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseEquationOfState::~multiphaseEquationOfState()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiphaseEquationOfState::updateP()
{
    volScalarField pByGamma
    (
        IOobject
        (
            "pRef",
            e_.mesh().time().timeName(),
            e_.mesh()
        ),
        e_.mesh(),
        dimensionedScalar("0", dimPressure, 0.0)
    );
    volScalarField rGamma
    (
        IOobject
        (
            "rGamma",
            e_.mesh().time().timeName(),
            e_.mesh()
        ),
        e_.mesh(),
        0.0
    );

    forAll(phases_, phasei)
    {
        eos_[phasei].update();
        volScalarField alphaByGamma
        (
            eos_[phasei]/(eos_[phasei].Gamma() - 1.0)
        );
        rGamma += alphaByGamma;
        pByGamma += eos_[phasei]*eos_[phasei].pByGamma();
    }

    p_ = pByGamma/rGamma;
    p_.max(small);
    p_.correctBoundaryConditions();
}

Foam::tmp<Foam::volScalarField>
Foam::multiphaseEquationOfState::c() const
{
    volScalarField alphaXiRhoCSqr
    (
        IOobject
        (
            "alphaXiRhoCSqr",
            e_.mesh().time().timeName(),
            e_.mesh()
        ),
        e_.mesh(),
        dimensionedScalar("0", dimPressure, 0.0)
    );
    volScalarField Xi
    (
        IOobject
        (
            "Xi",
            e_.mesh().time().timeName(),
            e_.mesh()
        ),
        e_.mesh(),
        0.0
    );

    forAll(eos_, phasei)
    {
        tmp<volScalarField> alphaXi(alphas_[phasei]*eos_[phasei].xi());
        tmp<volScalarField> h
        (
            (eos_[phasei].Gamma()*p_ + eos_[phasei].Pi())
           /((eos_[phasei].Gamma() - 1.0)*eos_[phasei].stabilizeRho())
        );
        tmp<volScalarField> cSqr
        (
            (h - eos_[phasei].delta())/eos_[phasei].xi()
        );

        Xi += alphaXi();
        alphaXiRhoCSqr += alphaXi*eos_[phasei].rho()*cSqr;
    }

    tmp<volScalarField> cSqr(alphaXiRhoCSqr/(rho_*Xi));
    cSqr.ref().max(small);
    return sqrt(cSqr);
}

// ************************************************************************* //

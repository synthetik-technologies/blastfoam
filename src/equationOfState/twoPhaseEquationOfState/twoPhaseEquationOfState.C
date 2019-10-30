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

#include "twoPhaseEquationOfState.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseEquationOfState::twoPhaseEquationOfState
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
    eos1_(equationOfState::New(e_, dict.subDict(phases_[0]))),
    alpha_(eos1_()),
    alphaRho1_(eos1_->alphaRho()),
    eos2_(equationOfState::New(e_, dict.subDict(phases_[1]))),
    alphaRho2_(eos2_->alphaRho())
{
    volScalarField& alpha2 = eos2_();
    alpha2 = 1.0 - alpha_;
    alphaRho2_ = alpha2*rho2();
    rho_ = alphaRho1_ + alphaRho2_;

    if (max(e_).value() < small)
    {
        e_ = alpha_*eos1_->e(p_) + (1.0 - alpha_)*eos2_->e(p_);
        e_.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseEquationOfState::~twoPhaseEquationOfState()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoPhaseEquationOfState::updateP()
{
    eos1_->update();
    eos2_->update();

    volScalarField rGamma
    (
        alpha_/(eos1_->Gamma() - 1.0) + (1.0 - alpha_)/(eos2_->Gamma() - 1.0)
    );

    p_ = (alpha_*eos1_->pByGamma() + (1.0 - alpha_)*eos2_->pByGamma())/rGamma;
    p_.max(small);
    p_.correctBoundaryConditions();
}

Foam::tmp<Foam::volScalarField>
Foam::twoPhaseEquationOfState::c() const
{
    volScalarField alphaXi1(alpha_*eos1_->xi());
    tmp<volScalarField> h1
    (
        (eos1_->Gamma()*p_ + eos1_->Pi())
       /((eos1_->Gamma() - 1.0)*eos1_->stabilizeRho())
    );
    tmp<volScalarField> c1Sqr((h1 - eos1_->delta())/eos1_->xi());

    volScalarField alphaXi2((1.0 - alpha_)*eos2_->xi());
    tmp<volScalarField> h2
    (
        (eos2_->Gamma()*p_ + eos2_->Pi())
       /((eos2_->Gamma() - 1.0)*eos2_->stabilizeRho())
    );
    tmp<volScalarField> c2Sqr((h2 - eos2_->delta())/eos2_->xi());

    volScalarField cSqr
    (
        (alphaXi1*rho1()*c1Sqr + alphaXi2*rho2()*c2Sqr)
        /(rho_*(alphaXi1 + alphaXi2))
    );
    cSqr.max(small);
    return sqrt(cSqr);
}

// ************************************************************************* //

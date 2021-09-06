/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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

#include "twoPhaseFluidBlastThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseFluidBlastThermo, 0);
    addToRunTimeSelectionTable
    (
        fluidBlastThermo,
        twoPhaseFluidBlastThermo,
        phase
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::twoPhaseBlastFluidMixture::HE
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (alpha1i_ < thermo1_.residualAlpha().value())
    {
        return mixture2Ptr_->HE(rho2i_, e, T);
    }
    else if (alpha2i_ < thermo2_.residualAlpha().value())
    {
        return mixture1Ptr_->HE(rho1i_, e, T);
    }
    return
        mixture1Ptr_->HE(rho1i_, e, T)*alpha1i_
      + mixture2Ptr_->HE(rho2i_, e, T)*alpha2i_;
}


Foam::scalar Foam::twoPhaseBlastFluidMixture::TRhoE
(
    const scalar T,
    const scalar rho,
    const scalar e
) const
{
    if (alpha1i_ < thermo1_.residualAlpha().value())
    {
        return mixture2Ptr_->TRhoE(T, rho2i_, e);
    }
    else if (alpha2i_ < thermo2_.residualAlpha().value())
    {
        return mixture1Ptr_->TRhoE(T, rho1i_, e);
    }
    return
        mixture1Ptr_->TRhoE(T, rho1i_, e)*alpha1i_
      + mixture2Ptr_->TRhoE(T, rho2i_, e)*alpha2i_;
}


Foam::scalar Foam::twoPhaseBlastFluidMixture::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (alpha1i_ < thermo1_.residualAlpha().value())
    {
        return mixture2Ptr_->Cp(rho2i_, e, T);
    }
    else if (alpha2i_ < thermo2_.residualAlpha().value())
    {
        return mixture1Ptr_->Cp(rho1i_, e, T);
    }
    return
        mixture1Ptr_->Cp(rho1i_, e, T)*alpha1i_
      + mixture2Ptr_->Cp(rho2i_, e, T)*alpha2i_;
}


Foam::scalar Foam::twoPhaseBlastFluidMixture::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (alpha1i_ < thermo1_.residualAlpha().value())
    {
        return mixture2Ptr_->Cv(rho2i_, e, T);
    }
    else if (alpha2i_ < thermo2_.residualAlpha().value())
    {
        return mixture1Ptr_->Cv(rho1i_, e, T);
    }
    return
        mixture1Ptr_->Cv(rho1i_, e, T)*alpha1i_
      + mixture2Ptr_->Cv(rho2i_, e, T)*alpha2i_;
}


Foam::scalar Foam::twoPhaseBlastFluidMixture::kappa
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (alpha1i_ < thermo1_.residualAlpha().value())
    {
        return mixture2Ptr_->kappa(rho2i_, e, T);
    }
    else if (alpha2i_ < thermo2_.residualAlpha().value())
    {
        return mixture1Ptr_->kappa(rho1i_, e, T);
    }
    return
        mixture1Ptr_->kappa(rho1i_, e, T)*alpha1i_
      + mixture2Ptr_->kappa(rho2i_, e, T)*alpha2i_;
}


Foam::scalar Foam::twoPhaseBlastFluidMixture::pRhoT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (alpha1i_ < thermo1_.residualAlpha().value())
    {
        return mixture2Ptr_->pRhoT(rho2i_, e, T);
    }
    else if (alpha2i_ < thermo2_.residualAlpha().value())
    {
        return mixture1Ptr_->pRhoT(rho1i_, e, T);
    }
    scalar alphaXi1(alpha1i_/(mixture1Ptr_->Gamma(rho1i_, e, T) - 1.0));
    scalar alphaXi2(alpha2i_/(mixture2Ptr_->Gamma(rho2i_, e, T) - 1.0));
    return
        (
            mixture1Ptr_->pRhoT(rho1i_, e, T)*alphaXi1
          + mixture2Ptr_->pRhoT(rho2i_, e, T)*alphaXi2
        )/(alphaXi1 + alphaXi2);
}


Foam::scalar Foam::twoPhaseBlastFluidMixture::mu
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (alpha1i_ < thermo1_.residualAlpha().value())
    {
        return mixture2Ptr_->mu(rho2i_, e, T);
    }
    else if (alpha2i_ < thermo2_.residualAlpha().value())
    {
        return mixture1Ptr_->mu(rho1i_, e, T);
    }
    return
        mixture1Ptr_->mu(rho1i_, e, T)*alpha1i_
      + mixture2Ptr_->mu(rho2i_, e, T)*alpha2i_;
}


Foam::scalar Foam::twoPhaseBlastFluidMixture::cSqr
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (alpha1i_ < thermo1_.residualAlpha().value())
    {
        return mixture2Ptr_->cSqr(p, rho2i_, e, T);
    }
    else if (alpha2i_ < thermo2_.residualAlpha().value())
    {
        return mixture1Ptr_->cSqr(p, rho1i_, e, T);
    }
    scalar alphaXi1(alpha1i_/(mixture1Ptr_->Gamma(rho1i_, e, T) - 1.0));
    scalar alphaXi2(alpha2i_/(mixture2Ptr_->Gamma(rho2i_, e, T) - 1.0));
    return
        (
            mixture1Ptr_->cSqr(p, rho1i_, e, T)*alphaXi1*rho1i_
          + mixture2Ptr_->cSqr(p, rho2i_, e, T)*alphaXi2*rho2i_
        )/(alphaXi1 + alphaXi2)/rho;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseFluidBlastThermo::twoPhaseFluidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    fluidBlastThermo(mesh, dict, phaseName),
    twoPhaseMixture(mesh, dict),
    rho1_
    (
        IOobject
        (
            IOobject::groupName("rho", phase1Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    thermo1_
    (
        phaseFluidBlastThermo::New
        (
            mesh,
            dict.subDict(phase1Name_),
            phase1Name_,
            phaseName
        )
    ),
    rho2_
    (
        IOobject
        (
            IOobject::groupName("rho", phase2Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    thermo2_
    (
        phaseFluidBlastThermo::New
        (
            mesh,
            dict.subDict(phase2Name_),
            phase2Name_,
            phaseName
        )
    ),
    mixture_(alpha1_, thermo1_(), thermo2_())
{
    //- Force reading of residual values
    thermo1_->read(dict.subDict(phase2Name_));
    thermo2_->read(dict.subDict(phase2Name_));

    this->residualAlpha_ =
        max(thermo1_->residualAlpha(), thermo2_->residualAlpha());
    this->residualRho_ = max(thermo1_->residualRho(), thermo2_->residualRho());

    // Update total density
    this->rho_ = alpha1_*rho1_ + alpha2_*rho2_;

    initializeFields();
}

void Foam::twoPhaseFluidBlastThermo::initializeModels()
{
    thermo1_->initializeModels();
    thermo2_->initializeModels();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseFluidBlastThermo::~twoPhaseFluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoPhaseFluidBlastThermo::read(const dictionary& dict)
{
    thermo1_->read(dict.subDict(thermo1_->name()));
    thermo2_->read(dict.subDict(thermo2_->name()));
}


void Foam::twoPhaseFluidBlastThermo::postUpdate()
{
    thermo1_->postUpdate();
    thermo2_->postUpdate();
}


void Foam::twoPhaseFluidBlastThermo::solve()
{
    thermo1_->solve();
    thermo2_->solve();
}


void Foam::twoPhaseFluidBlastThermo::update()
{

    thermo1_->update();
    thermo2_->update();
}


void Foam::twoPhaseFluidBlastThermo::updateRho(const volScalarField& p)
{
    thermo1_->updateRho(p);
    thermo2_->updateRho(p);
    this->rho_ = alpha1_*thermo1_->rho() + alpha2_*thermo2_->rho();
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseFluidBlastThermo::ESource() const
{
    return alpha1_*thermo1_->ESource() + alpha2_*thermo2_->ESource();
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseFluidBlastThermo::calce(const volScalarField& p) const
{
    return alpha1_*thermo1_->calce(p) + alpha2_*thermo2_->calce(p);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellpRhoT(const label celli) const
{
    scalar alphaXi1
    (
        alpha1_[celli]/(thermo1_->cellGamma(celli) - 1.0)
    );
    scalar alphaXi2
    (
        alpha2_[celli]/(thermo2_->cellGamma(celli) - 1.0)
    );

    return
        (
            alphaXi1*thermo1_->cellpRhoT(celli)*pos(alphaXi1 - small)
          + alphaXi2*thermo2_->cellpRhoT(celli)*pos(alphaXi2 - small)
        )/(alphaXi1 + alphaXi2);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::celldpdRho(const label celli) const
{
    scalar alphaXi1
    (
        alpha1_[celli]/(thermo1_->cellGamma(celli) - 1.0)
    );
    scalar alphaXi2
    (
        alpha2_[celli]/(thermo2_->cellGamma(celli) - 1.0)
    );

    return
        (
            alphaXi1*thermo1_->celldpdRho(celli)*pos(alphaXi1 - small)
          + alphaXi2*thermo2_->celldpdRho(celli)*pos(alphaXi2 - small)
        )/(alphaXi1 + alphaXi2);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::celldpde(const label celli) const
{
    scalar alphaXi1
    (
        alpha1_[celli]/(thermo1_->cellGamma(celli) - 1.0)
    );
    scalar alphaXi2
    (
        alpha2_[celli]/(thermo2_->cellGamma(celli) - 1.0)
    );

    return
        (
            alphaXi1*thermo1_->celldpde(celli)*pos(alphaXi1 - small)
          + alphaXi2*thermo2_->celldpde(celli)*pos(alphaXi2 - small)
        )/(alphaXi1 + alphaXi2);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellGamma(const label celli) const
{
    return
        alpha1_[celli]*thermo1_->cellGamma(celli)
      + alpha2_[celli]*thermo2_->cellGamma(celli);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return
        alpha1_*thermo1_->he(p, T)
      + alpha2_*thermo2_->he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    return
        UIndirectList<scalar>(alpha1_(), cells)()*thermo1_->he(T, cells)
      + UIndirectList<scalar>(alpha2_(), cells)()*thermo2_->he(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::he
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->he(T, patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->he(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::hs() const
{
    return alpha1_*thermo1_->hs() + alpha2_*thermo2_->hs();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return alpha1_*thermo1_->hs(p, T) + alpha2_*thermo2_->hs(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    return
        UIndirectList<scalar>(alpha1_(), cells)()*thermo1_->hs(T, cells)
      + UIndirectList<scalar>(alpha2_(), cells)()*thermo2_->hs(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->hs(T, patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->hs(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::ha() const
{
    return alpha1_*thermo1_->ha() + alpha2_*thermo2_->ha();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return alpha1_*thermo1_->ha(p, T) + alpha2_*thermo2_->ha(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    return
        UIndirectList<scalar>(alpha1_(), cells)()*thermo1_->ha(T, cells)
      + UIndirectList<scalar>(alpha2_(), cells)()*thermo2_->ha(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::ha
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->ha(T, patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->ha(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::hc() const
{
    return
        alpha1_*thermo1_->hc()
      + alpha2_*thermo2_->hc();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::flameT() const
{
    return
        alpha1_*thermo1_->flameT()
      + alpha2_*thermo2_->flameT();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::THE() const
{
    return
        alpha1_*thermo1_->THE()
      + alpha2_*thermo2_->THE();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::THE
(
    const volScalarField& he,
    const volScalarField& p,
    const volScalarField& T0
) const
{
    return
        alpha1_*thermo1_->THE(he, p, T0)
      + alpha2_*thermo2_->THE(he, p, T0);
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidBlastThermo::THE
(
    const scalarField& he,
    const scalarField& T,
    const labelList& cells
) const
{
    return
        UIndirectList<scalar>(alpha1_(), cells)()*thermo1_->THE(he, T, cells)
      + UIndirectList<scalar>(alpha2_(), cells)()*thermo2_->THE(he, T, cells);
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidBlastThermo::THE
(
    const scalarField& he,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->THE(he, T, patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->THE(he, T, patchi);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellTHE
(
    const scalar he,
    const scalar T,
    const label celli
) const
{
    return
        alpha1_[celli]*thermo1_->cellTHE(he, T, celli)
      + alpha2_[celli]*thermo2_->cellTHE(he, T, celli);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->Cp(T, patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->Cp(T, patchi);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellCp
(
    const scalar T,
    const label celli
) const
{
    return
        alpha1_[celli]*thermo1_->cellCp(T, celli)
      + alpha2_[celli]*thermo2_->cellCp(T, celli);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->Cv(T, patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->Cv(T, patchi);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellCv
(
    const scalar T,
    const label celli
) const
{
    return
        alpha1_[celli]*thermo1_->cellCv(T, celli)
      + alpha2_[celli]*thermo2_->cellCv(T, celli);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->Cpv(T, patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->Cpv(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::W() const
{
    return alpha1_*thermo1_->W() + alpha2_*thermo2_->W();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::W
(
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->W(patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->W(patchi);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellW(const label celli) const
{
    return
        alpha1_[celli]*thermo1_->cellW(celli)
      + alpha2_[celli]*thermo2_->cellW(celli);
}

// ************************************************************************* //

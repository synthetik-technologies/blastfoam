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

void Foam::twoPhaseFluidBlastThermo::calculate()
{
    volScalarField T0(T_);
    T_ = Zero;
    thermo1_->calculateTemperature(alpha1_, this->he(), T0, T_);
    thermo2_->calculateTemperature(alpha2_, this->he(), T0, T_);
    T_.correctBoundaryConditions();

    if (min(T_).value() < this->TLow_)
    {
        T0 = T_;
        T_.max(this->TLow_);
        volScalarField he0(this->he());
        this->he() = Zero;
        thermo1_->calculateEnergy(alpha1_, T0, he0, T_, this->he());
        thermo2_->calculateEnergy(alpha2_, T0, he0, T_, this->he());
    }
    this->he().correctBoundaryConditions();

    volScalarField XiSum
    (
        volScalarField::New("XiSum", mesh(), dimensionedScalar(dimless, 0.0))
    );
    volScalarField pXiSum
    (
        volScalarField::New
        (
            "pXiSum",
            mesh(),
            dimensionedScalar(dimPressure, 0.0)
        )
    );
    volScalarField cSqrRhoXiSum
    (
        volScalarField::New
        (
            "cSqrRhoXiSum",
            mesh(),
            dimensionedScalar(sqr(dimVelocity)*dimDensity, 0.0)
        )
    );

    this->Cp_ = Zero;
    this->Cv_ = Zero;
    this->alpha_ = Zero;
    this->mu_ = Zero;

    thermo1_->calculate
    (
        alpha1_,
        this->he(),
        this->T_,
        this->Cp_,
        this->Cv_,
        this->alpha_,
        this->mu_,
        pXiSum,
        XiSum
    );
    thermo2_->calculate
    (
        alpha2_,
        this->he(),
        this->T_,
        this->Cp_,
        this->Cv_,
        this->alpha_,
        this->mu_,
        pXiSum,
        XiSum
    );
    XiSum.max(small);
    this->p_ = pXiSum/XiSum;
    this->p_.correctBoundaryConditions();

    thermo1_->calculateSpeedOfSound
    (
        alpha1_,
        cSqrRhoXiSum
    );
    thermo2_->calculateSpeedOfSound
    (
        alpha2_,
        cSqrRhoXiSum
    );

    cSqrRhoXiSum.max(small);
    this->speedOfSound_ =
        sqrt
        (
            cSqrRhoXiSum
           /max(this->rho_*XiSum, this->residualRho_)
        );
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
    )
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


void Foam::twoPhaseFluidBlastThermo::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
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

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
#include "equation.H"
#include "NewtonRaphsonUnivariateRootSolver.H"
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

// * * * * * * * * * * twoPhaseEEquation Member Functions  * * * * * * * * * //

namespace Foam
{
    class twoPhaseEEquation
    :
        public equation
    {
        const twoPhaseFluidBlastThermo& thermo_;
        const volScalarField& p_;
        mutable volScalarField* e_;

    public:
        twoPhaseEEquation(twoPhaseFluidBlastThermo& thermo)
        :
            equation(-great, great),
            thermo_(thermo),
            p_(thermo_.p()),
            e_(&thermo.he())
        {}
        virtual label nDerivatives() const
        {
            return 1;
        }
        virtual scalar fx(const scalar e, const label li) const
        {
            (*e_)[li] = e;
            return thermo_.cellpRhoT(li, false) - p_[li];
        }
        virtual scalar dfdx(const scalar e, const label li) const
        {
            (*e_)[li] = e;
            return thermo_.celldpde(li);
        }
    };


    class twoPhaseTHEEquation
    :
        public equation
    {
        const twoPhaseFluidBlastThermo& thermo_;
        const volScalarField& he_;
        label patchi_;

    public:
        twoPhaseTHEEquation
        (
            twoPhaseFluidBlastThermo& thermo,
            const scalar TLow
        )
        :
            equation(-great, great),
            thermo_(thermo),
            he_(thermo.he()),
            patchi_(-1)
        {}

        label& patch()
        {
            return patchi_;
        }

        virtual label nDerivatives() const
        {
            return 1;
        }
        virtual scalar fx(const scalar T, const label li) const
        {
            return
                patchi_ == -1
              ? thermo_.cellHE(T, li) - he_[li]
              : thermo_.patchFaceHE(T, patchi_, li)
              - he_.boundaryField()[patchi_][li];
        }
        virtual scalar dfdx(const scalar T, const label li) const
        {
            return
                patchi_ == -1
              ? thermo_.cellCpv(T, li)
              : thermo_.patchFaceCpv(T, patchi_, li);
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::twoPhaseFluidBlastThermo::calculate()
{
    scalarField& TCells = T_.primitiveFieldRef();
    scalarField& heCells = this->he().primitiveFieldRef();
    volScalarField::Boundary& bT = T_.boundaryFieldRef();
    volScalarField::Boundary& bhe = this->he().boundaryFieldRef();

    twoPhaseTHEEquation eqn(*this, this->TLow_);
    NewtonRaphsonUnivariateRootSolver solver(eqn, dictionary());
    forAll(TCells, celli)
    {
        TCells[celli] = solver.solve(TCells[celli], celli);
    }
    forAll(bT, patchi)
    {
        eqn.patch() = patchi;
        scalarField& pT = bT[patchi];
        forAll(pT, facei)
        {
            pT[facei] = solver.solve(pT[facei], facei);
        }
    }
    T_.correctBoundaryConditions();

    if (min(T_).value() <= this->TLow_)
    {
        forAll(TCells, celli)
        {
            if (TCells[celli] <= this->TLow_)
            {
                TCells[celli] = this->TLow_;
                heCells[celli] = this->cellHE(this->TLow_, celli);
            }
        }
        forAll(bhe, patchi)
        {
            scalarField& pT = bT[patchi];
            scalarField& phe = bhe[patchi];

            pT = max(pT, this->TLow_);
            phe = this->he(pT, patchi);
        }
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

bool Foam::twoPhaseFluidBlastThermo::read()
{
    thermo1_->read(this->subDict(thermo1_->name()));
    thermo2_->read(this->subDict(thermo2_->name()));
    return true;
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
    tmp<volScalarField> eInitTmp(volScalarField::New("eInit", e_));
    volScalarField& eInit(eInitTmp.ref());
    twoPhaseFluidBlastThermo& thermo
    (
        const_cast<twoPhaseFluidBlastThermo&>(*this)
    );
    volScalarField& e(const_cast<volScalarField&>(e_));

    twoPhaseEEquation eqn(thermo);
    dictionary dict;
    dict.add("tolerance", 1e-6);
    NewtonRaphsonUnivariateRootSolver solver(eqn, dict);
    forAll(eInit, celli)
    {
        if (mag(celldpde(celli)) < small)
        {
            eInit[celli] = cellHE(T_[celli], celli);
        }
        else
        {
            const scalar e0 = e[celli];
            eInit[celli] = solver.solve(e0, celli);
            e[celli] = e0;
        }
    }
    eInit +=
        alpha1_*thermo1_->initESource() + alpha2_*thermo2_->initESource();

    return eInitTmp;
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellpRhoT
(
    const label celli,
    const bool limit
) const
{
    if (alpha2_[celli] < thermo2_->residualAlpha().value())
    {
        return thermo1_->cellpRhoT(celli, limit);
    }
    if (alpha1_[celli] < thermo1_->residualAlpha().value())
    {
        return thermo2_->cellpRhoT(celli, limit);
    }
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
            alphaXi1*thermo1_->cellpRhoT(celli, limit)
          + alphaXi2*thermo2_->cellpRhoT(celli, limit)
        )/(alphaXi1 + alphaXi2);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::celldpdRho(const label celli) const
{
    if (alpha2_[celli] < thermo2_->residualAlpha().value())
    {
        return thermo1_->celldpdRho(celli);
    }
    if (alpha1_[celli] < thermo1_->residualAlpha().value())
    {
        return thermo2_->celldpdRho(celli);
    }
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
            alphaXi1*thermo1_->celldpdRho(celli)
          + alphaXi2*thermo2_->celldpdRho(celli)
        )/(alphaXi1 + alphaXi2);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::celldpde(const label celli) const
{
    if (alpha2_[celli] < thermo2_->residualAlpha().value())
    {
        return thermo1_->celldpde(celli);
    }
    if (alpha1_[celli] < thermo1_->residualAlpha().value())
    {
        return thermo2_->celldpde(celli);
    }
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
            alphaXi1*thermo1_->celldpde(celli)
          + alphaXi2*thermo2_->celldpde(celli)
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


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellHE
(
    const scalar T,
    const label celli
) const
{
    if (alpha2_[celli] < residualAlpha_.value())
    {
        return thermo1_->cellHE(T, celli);
    }
    else if (alpha1_[celli] < residualAlpha_.value())
    {
        return thermo2_->cellHE(T, celli);
    }
    scalar alphaRho1 = alpha1_[celli]*rho1_[celli];
    scalar alphaRho2 = alpha2_[celli]*rho2_[celli];
    return
        (
            alphaRho1*thermo1_->cellHE(T, celli)
          + alphaRho2*thermo2_->cellHE(T, celli)
        )/(alphaRho1 + alphaRho2);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::patchFaceHE
(
    const scalar T,
    const label patchi,
    const label facei
) const
{
    if (alpha2_.boundaryField()[patchi][facei] < residualAlpha_.value())
    {
        return thermo1_->patchFaceHE(T, patchi, facei);
    }
    if (alpha1_.boundaryField()[patchi][facei] < residualAlpha_.value())
    {
        return thermo2_->patchFaceHE(T, patchi, facei);
    }
    scalar alphaRho1 =
        alpha1_.boundaryField()[patchi][facei]
       *rho1_.boundaryField()[patchi][facei];
    scalar alphaRho2 =
        alpha2_.boundaryField()[patchi][facei]
       *rho2_.boundaryField()[patchi][facei];
    return
        (
            alphaRho1*thermo1_->patchFaceHE(T, patchi, facei)
          + alphaRho2*thermo2_->patchFaceHE(T, patchi, facei)
        )/(alphaRho1 + alphaRho2);
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


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellCpv
(
    const scalar T,
    const label celli
) const
{
    if (alpha2_[celli] < residualAlpha_.value())
    {
        return thermo1_->cellCpv(T, celli);
    }
    else if (alpha1_[celli] < residualAlpha_.value())
    {
        return thermo2_->cellCpv(T, celli);
    }
    scalar alphaRho1 = alpha1_[celli]*rho1_[celli];
    scalar alphaRho2 = alpha2_[celli]*rho2_[celli];
    return
        (
            alphaRho1*thermo1_->cellCpv(T, celli)
          + alphaRho2*thermo2_->cellCpv(T, celli)
        )/(alphaRho1 + alphaRho2);
}



Foam::scalar Foam::twoPhaseFluidBlastThermo::patchFaceCpv
(
    const scalar T,
    const label patchi,
    const label facei
) const
{
    if (alpha2_.boundaryField()[patchi][facei] < residualAlpha_.value())
    {
        return thermo1_->patchFaceCpv(T, patchi, facei);
    }
    if (alpha1_.boundaryField()[patchi][facei] < residualAlpha_.value())
    {
        return thermo2_->patchFaceCpv(T, patchi, facei);
    }
    scalar alphaRho1 =
        alpha1_.boundaryField()[patchi][facei]
       *rho1_.boundaryField()[patchi][facei];
    scalar alphaRho2 =
        alpha2_.boundaryField()[patchi][facei]
       *rho2_.boundaryField()[patchi][facei];

    return
        (
            alphaRho1*thermo1_->patchFaceCpv(T, patchi, facei)
          + alphaRho2*thermo2_->patchFaceCpv(T, patchi, facei)
        )/(alphaRho1 + alphaRho2);
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

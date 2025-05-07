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
#include "Equation.H"
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::twoPhaseFluidBlastThermo::calculate()
{
    scalarField& TCells = T_.primitiveFieldRef();
    scalarField& heCells = this->he().primitiveFieldRef();

    THEEqn_.setCells();
    forAll(TCells, celli)
    {
        TCells[celli] = THESolver_->solve(TCells[celli], celli);
        if (TCells[celli] <= this->TLow_)
        {
            TCells[celli] = this->TLow_;
            heCells[celli] = this->cellhe(this->TLow_, celli);
        }
    }


    volScalarField::Boundary& bT = T_.boundaryFieldRef();
    volScalarField::Boundary& bhe = this->he().boundaryFieldRef();
    forAll(bT, patchi)
    {
        THEEqn_.setPatch(patchi);
        fvPatchScalarField& pT = bT[patchi];
        fvPatchScalarField& phe = bhe[patchi];
        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                phe[facei] = this->patchFacehe
                (
                    pT[facei],
                    patchi,
                    facei
                );
            }
        }
        else
        {
            forAll(pT, facei)
            {
                pT[facei] = THESolver_->solve(pT[facei], facei);
                if (pT[facei] < TLow_)
                {
                    pT[facei] = TLow_;
                    phe[facei] = this->patchFacehe
                    (
                        pT[facei],
                        patchi,
                        facei
                    );
                }
            }
        }
    }

    volScalarField XiSum
    (
        volScalarField::New
        (
            "XiSum",
            twoPhases::mesh(),
            dimensionedScalar(dimless, 0.0)
        )
    );
    volScalarField pXiSum
    (
        volScalarField::New
        (
            "pXiSum",
            twoPhases::mesh(),
            dimensionedScalar(dimPressure, 0.0)
        )
    );
    volScalarField cSqrRhoXiSum
    (
        volScalarField::New
        (
            "cSqrRhoXiSum",
            twoPhases::mesh(),
            dimensionedScalar(sqr(dimVelocity)*dimDensity, 0.0)
        )
    );

    this->Cp_ = Zero;
    this->Cv_ = Zero;
    this->kappa_ = Zero;
    this->mu_ = Zero;

    thermo1_->calculate
    (
        this->alpha1(),
        this->he(),
        this->T_,
        this->Cp_,
        this->Cv_,
        this->kappa_,
        this->mu_,
        pXiSum,
        XiSum
    );
    thermo2_->calculate
    (
        this->alpha2(),
        this->he(),
        this->T_,
        this->Cp_,
        this->Cv_,
        this->kappa_,
        this->mu_,
        pXiSum,
        XiSum
    );
    XiSum.max(small);
    this->p_ = pXiSum/XiSum;
    this->p_.correctBoundaryConditions();

    thermo1_->calculateSpeedOfSound
    (
        this->alpha1(),
        cSqrRhoXiSum
    );
    thermo2_->calculateSpeedOfSound
    (
        this->alpha2(),
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
    fluidBlastThermo(mesh, dict, phaseName, word::null, false),
    phase1Name_(wordList(dict.lookup("phases"))[0]),
    alpha1_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase1Name_),
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    rho1_
    (
        IOobject
        (
            IOobject::groupName("rho", phase1Name_),
            mesh.time().name(),
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
    phase2Name_(wordList(dict.lookup("phases"))[1]),
    alpha2_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase2Name_),
            mesh.time().name(),
            mesh
        ),
        1.0 - alpha1_
    ),
    rho2_
    (
        IOobject
        (
            IOobject::groupName("rho", phase2Name_),
            mesh.time().name(),
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
    TEqn_(*this, dict.optionalSubDict("eSolverCoeffs")),
    TSolver_(nullptr),
    THEEqn_(*this, this->TLow_, dict.optionalSubDict("TSolverCoeffs")),
    THESolver_(nullptr)
{
    // Select the solvers for energy and temperature
    if (dict.isDict("eSolverCoeffs"))
    {
        const dictionary& eDict(dict.subDict("eSolverCoeffs"));
        TSolver_ =
            univariateRootSolver::New
            (
                eDict.lookupOrDefault
                (
                    "solver",
                    NewtonRaphsonUnivariateRootSolver::typeName
                ),
                TEqn_,
                eDict
            );
    }
    else
    {
        TSolver_ =
            univariateRootSolver::New
            (
                NewtonRaphsonUnivariateRootSolver::typeName,
                TEqn_,
                dict
            );
    }
    if (dict.isDict("TSolverCoeffs"))
    {
        const dictionary& TDict(dict.subDict("TSolverCoeffs"));
        THESolver_ =
            univariateRootSolver::New
            (
                TDict.lookupOrDefault
                (
                    "solver",
                    NewtonRaphsonUnivariateRootSolver::typeName
                ),
                THEEqn_,
                TDict
            );
    }
    else
    {
        THESolver_ =
            univariateRootSolver::New
            (
                NewtonRaphsonUnivariateRootSolver::typeName,
                THEEqn_,
                dict
            );
    }

    //- Force reading of residual values
    thermo1_->read(dict.subDict(phase2Name_));
    thermo2_->read(dict.subDict(phase2Name_));

    this->residualAlpha_ =
        max(thermo1_->residualAlpha(), thermo2_->residualAlpha());
    this->residualRho_ = max(thermo1_->residualRho(), thermo2_->residualRho());

    // Update total density
    this->rho_ = this->alpha1()*rho1_ + this->alpha2()*rho2_;

    // Initial guess for e
    if (!this->e_.headerOk())
    {
        this->e_ == he(p_, T_);
    }

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
    thermo1_->updateRho(this->alpha1(), p);
    thermo2_->updateRho(this->alpha2(), p);
    this->rho_ = this->alpha1()*thermo1_->rho() + this->alpha2()*thermo2_->rho();
    this->rho_.correctBoundaryConditions();
}


bool Foam::twoPhaseFluidBlastThermo::contains(const word& specieName) const
{
    return thermo1_->contains(specieName) || thermo2_->contains(specieName);
}

void Foam::twoPhaseFluidBlastThermo::addDelta
(
    const word& name,
    tmp<volScalarField>& delta
)
{
    if (thermo1_->contains(name))
    {
        thermo1_->addDelta(name, delta);
    }
    else if (thermo2_->contains(name))
    {
        thermo2_->addDelta(name, delta);
    }
}


void Foam::twoPhaseFluidBlastThermo::addDelta
(
    const word& name,
    const volScalarField::Internal& delta
)
{
    if (thermo1_->contains(name))
    {
        thermo1_->addDelta(name, delta);
    }
    else if (thermo2_->contains(name))
    {
        thermo2_->addDelta(name, delta);
    }
}


void Foam::twoPhaseFluidBlastThermo::addSource
(
    const word& name,
    tmp<fvScalarMatrix>& source
)
{
    if (thermo1_->contains(name))
    {
        thermo1_->addSource(name, source);
    }
    else if (thermo2_->contains(name))
    {
        thermo2_->addSource(name, source);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseFluidBlastThermo::ESource() const
{
    return this->alpha1()*thermo1_->ESource() + this->alpha2()*thermo2_->ESource();
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseFluidBlastThermo::calce(const volScalarField& p) const
{
    tmp<volScalarField> eInitTmp(volScalarField::New("eInit", e_));
    volScalarField& eInit(eInitTmp.ref());

    if (!this->rho_.time().restart())
    {
        forAll(eInit, celli)
        {
            scalar Tinit = this->T_[celli];
            if (mag(celldpdT(celli)) > small)
            {
                TEqn_.save(p[celli], celli);
                Tinit = TSolver_->solve(T_[celli], celli);
                TEqn_.reset(celli);
            }
            eInit[celli] = cellhe(Tinit, celli);
        }
        eInit +=
            this->alpha1()*thermo1_->initESource() + this->alpha2()*thermo2_->initESource();
    }
    else
    {
        forAll(eInit, celli)
        {
            eInit[celli] = cellhe(this->T_[celli], celli);
        }
    }

    return eInitTmp;
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::calcCelle
(
    const scalar p,
    const label celli
) const
{
    scalar Tinit = this->T_[celli];
    if (mag(celldpdT(celli)) > small)
    {
        TEqn_.save(p, celli);
        Tinit = TSolver_->solve(T_[celli], celli);
        TEqn_.reset(celli);
    }
    return cellhe(Tinit, celli);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellpRhoT
(
    const label celli,
    const bool limit
) const
{
    if (this->alpha2()[celli] < thermo2_->residualAlpha().value())
    {
        return thermo1_->cellpRhoT(celli, limit);
    }
    if (this->alpha1()[celli] < thermo1_->residualAlpha().value())
    {
        return thermo2_->cellpRhoT(celli, limit);
    }
    scalar alphaXi1
    (
        this->alpha1()[celli]/(thermo1_->cellGamma(celli) - 1.0)
    );
    scalar alphaXi2
    (
        this->alpha2()[celli]/(thermo2_->cellGamma(celli) - 1.0)
    );

    return
        (
            alphaXi1*thermo1_->cellpRhoT(celli, limit)
          + alphaXi2*thermo2_->cellpRhoT(celli, limit)
        )/(alphaXi1 + alphaXi2);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::patchFacepRhoT
(
    const label patchi,
    const label facei,
    const bool limit
) const
{
    const scalar a1 = alpha1_.boundaryField()[patchi][facei];
    const scalar a2 = alpha2_.boundaryField()[patchi][facei];
    if (a2 < thermo2_->residualAlpha().value())
    {
        return thermo1_->patchFacepRhoT(patchi, facei, limit);
    }
    if (a1 < thermo1_->residualAlpha().value())
    {
        return thermo2_->patchFacepRhoT(patchi, facei, limit);
    }
    const scalar alphaXi1 = a1/(thermo1_->patchFaceGamma(patchi, facei) - 1.0);
    const scalar alphaXi2 = a2/(thermo2_->patchFaceGamma(patchi, facei) - 1.0);

    return
        (
            alphaXi1*thermo1_->patchFacepRhoT(patchi, facei, limit)
          + alphaXi2*thermo2_->patchFacepRhoT(patchi, facei, limit)
        )/(alphaXi1 + alphaXi2);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::celldpdRho(const label celli) const
{
    if (this->alpha2()[celli] < thermo2_->residualAlpha().value())
    {
        return thermo1_->celldpdRho(celli);
    }
    if (this->alpha1()[celli] < thermo1_->residualAlpha().value())
    {
        return thermo2_->celldpdRho(celli);
    }
    scalar alphaXi1
    (
        this->alpha1()[celli]/(thermo1_->cellGamma(celli) - 1.0)
    );
    scalar alphaXi2
    (
        this->alpha2()[celli]/(thermo2_->cellGamma(celli) - 1.0)
    );

    return
        (
            alphaXi1*thermo1_->celldpdRho(celli)
          + alphaXi2*thermo2_->celldpdRho(celli)
        )/(alphaXi1 + alphaXi2);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::celldpde(const label celli) const
{
    if (this->alpha2()[celli] < thermo2_->residualAlpha().value())
    {
        return thermo1_->celldpde(celli);
    }
    if (this->alpha1()[celli] < thermo1_->residualAlpha().value())
    {
        return thermo2_->celldpde(celli);
    }
    scalar alphaXi1
    (
        this->alpha1()[celli]/(thermo1_->cellGamma(celli) - 1.0)
    );
    scalar alphaXi2
    (
        this->alpha2()[celli]/(thermo2_->cellGamma(celli) - 1.0)
    );

    return
        (
            alphaXi1*thermo1_->celldpde(celli)
          + alphaXi2*thermo2_->celldpde(celli)
        )/(alphaXi1 + alphaXi2);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::celldpdT(const label celli) const
{
    if (this->alpha2()[celli] < thermo2_->residualAlpha().value())
    {
        return thermo1_->celldpdT(celli);
    }
    if (this->alpha1()[celli] < thermo1_->residualAlpha().value())
    {
        return thermo2_->celldpdT(celli);
    }
    scalar alphaXi1
    (
        this->alpha1()[celli]/(thermo1_->cellGamma(celli) - 1.0)
    );
    scalar alphaXi2
    (
        this->alpha2()[celli]/(thermo2_->cellGamma(celli) - 1.0)
    );

    return
        (
            alphaXi1*thermo1_->celldpdT(celli)
          + alphaXi2*thermo2_->celldpdT(celli)
        )/(alphaXi1 + alphaXi2);
}

Foam::scalar Foam::twoPhaseFluidBlastThermo::cellGamma(const label celli) const
{
    return
        this->alpha1()[celli]*thermo1_->cellGamma(celli)
      + this->alpha2()[celli]*thermo2_->cellGamma(celli);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::patchFaceGamma
(
    const label patchi,
    const label facei
) const
{
    return
        alpha1_.boundaryField()[patchi][facei]*thermo1_->patchFaceGamma(patchi, facei)
      + alpha2_.boundaryField()[patchi][facei]*thermo2_->patchFaceGamma(patchi, facei);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return
        this->alpha1()*thermo1_->he(p, T)
      + this->alpha2()*thermo2_->he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    return
        UIndirectList<scalar>(this->alpha1()(), cells)()*thermo1_->he(T, cells)
      + UIndirectList<scalar>(this->alpha2()(), cells)()*thermo2_->he(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::he
(
    const scalarField& T,
    const label patchi
) const
{
    return
        this->alpha1().boundaryField()[patchi]*thermo1_->he(T, patchi)
      + this->alpha2().boundaryField()[patchi]*thermo2_->he(T, patchi);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::he
(
    const scalarField& T,
    const fvSource& source
) const
{
    return
        scalarField(this->alpha1(), source.cells())
       *thermo1_->he(T, source)
      + scalarField(this->alpha2(), source.cells())
       *thermo2_->he(T, source);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellhe
(
    const scalar T,
    const label celli
) const
{
    if (this->alpha2()[celli] < residualAlpha_.value())
    {
        return thermo1_->cellhe(T, celli);
    }
    else if (this->alpha1()[celli] < residualAlpha_.value())
    {
        return thermo2_->cellhe(T, celli);
    }
    scalar alphaRho1 = this->alpha1()[celli]*rho1_[celli];
    scalar alphaRho2 = this->alpha2()[celli]*rho2_[celli];
    return
        (
            alphaRho1*thermo1_->cellhe(T, celli)
          + alphaRho2*thermo2_->cellhe(T, celli)
        )/(alphaRho1 + alphaRho2);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::patchFacehe
(
    const scalar T,
    const label patchi,
    const label facei
) const
{
    if (this->alpha2().boundaryField()[patchi][facei] < residualAlpha_.value())
    {
        return thermo1_->patchFacehe(T, patchi, facei);
    }
    if (this->alpha1().boundaryField()[patchi][facei] < residualAlpha_.value())
    {
        return thermo2_->patchFacehe(T, patchi, facei);
    }
    scalar alphaRho1 =
        this->alpha1().boundaryField()[patchi][facei]
       *rho1_.boundaryField()[patchi][facei];
    scalar alphaRho2 =
        this->alpha2().boundaryField()[patchi][facei]
       *rho2_.boundaryField()[patchi][facei];
    return
        (
            alphaRho1*thermo1_->patchFacehe(T, patchi, facei)
          + alphaRho2*thermo2_->patchFacehe(T, patchi, facei)
        )/(alphaRho1 + alphaRho2);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::hs() const
{
    return this->alpha1()*thermo1_->hs() + this->alpha2()*thermo2_->hs();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return this->alpha1()*thermo1_->hs(p, T) + this->alpha2()*thermo2_->hs(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    return
        UIndirectList<scalar>(this->alpha1()(), cells)()*thermo1_->hs(T, cells)
      + UIndirectList<scalar>(this->alpha2()(), cells)()*thermo2_->hs(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    return
        this->alpha1().boundaryField()[patchi]*thermo1_->hs(T, patchi)
      + this->alpha2().boundaryField()[patchi]*thermo2_->hs(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::ha() const
{
    return this->alpha1()*thermo1_->ha() + this->alpha2()*thermo2_->ha();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return this->alpha1()*thermo1_->ha(p, T) + this->alpha2()*thermo2_->ha(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    return
        UIndirectList<scalar>(this->alpha1()(), cells)()*thermo1_->ha(T, cells)
      + UIndirectList<scalar>(this->alpha2()(), cells)()*thermo2_->ha(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::ha
(
    const scalarField& T,
    const label patchi
) const
{
    return
        this->alpha1().boundaryField()[patchi]*thermo1_->ha(T, patchi)
      + this->alpha2().boundaryField()[patchi]*thermo2_->ha(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::hc() const
{
    return
        this->alpha1()*thermo1_->hc()
      + this->alpha2()*thermo2_->hc();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::The() const
{
    return
        this->alpha1()*thermo1_->The()
      + this->alpha2()*thermo2_->The();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::The
(
    const volScalarField& he,
    const volScalarField& p,
    const volScalarField& T0
) const
{
    return
        this->alpha1()*thermo1_->The(he, p, T0)
      + this->alpha2()*thermo2_->The(he, p, T0);
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidBlastThermo::The
(
    const scalarField& he,
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(this->alpha1(), cells)*thermo1_->The(he, T, cells)
      + scalarField(this->alpha2(), cells)*thermo2_->The(he, T, cells);
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidBlastThermo::The
(
    const scalarField& he,
    const scalarField& T,
    const label patchi
) const
{
    return
        this->alpha1().boundaryField()[patchi]*thermo1_->The(he, T, patchi)
      + this->alpha2().boundaryField()[patchi]*thermo2_->The(he, T, patchi);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellThe
(
    const scalar he,
    const scalar T,
    const label celli
) const
{
    return
        this->alpha1()[celli]*thermo1_->cellThe(he, T, celli)
      + this->alpha2()[celli]*thermo2_->cellThe(he, T, celli);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return
        (
            this->alpha1().boundaryField()[patchi]
           *rho1_.boundaryField()[patchi]
           *thermo1_->Cp(T, patchi)
          + this->alpha2().boundaryField()[patchi]
           *rho2_.boundaryField()[patchi]
           *thermo2_->Cp(T, patchi)
        )/rho_.boundaryField()[patchi];
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellCp
(
    const scalar T,
    const label celli
) const
{
    return
        (
            this->alpha1()[celli]*rho1_[celli]*thermo1_->cellCp(T, celli)
          + this->alpha2()[celli]*rho2_[celli]*thermo2_->cellCp(T, celli)
        )/rho_[celli];
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return
        (
            this->alpha1().boundaryField()[patchi]
           *rho1_.boundaryField()[patchi]
           *thermo1_->Cv(T, patchi)
          + this->alpha2().boundaryField()[patchi]
           *rho2_.boundaryField()[patchi]
           *thermo2_->Cv(T, patchi)
        )/rho_.boundaryField()[patchi];
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellCv
(
    const scalar T,
    const label celli
) const
{
    return
        (
            this->alpha1()[celli]*rho1_[celli]*thermo1_->cellCv(T, celli)
          + this->alpha2()[celli]*rho2_[celli]*thermo2_->cellCv(T, celli)
        )/rho_[celli];
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return
        (
            this->alpha1().boundaryField()[patchi]
           *rho1_.boundaryField()[patchi]
           *thermo1_->Cpv(T, patchi)
          + this->alpha2().boundaryField()[patchi]
           *rho2_.boundaryField()[patchi]
           *thermo2_->Cpv(T, patchi)
        )/rho_.boundaryField()[patchi];
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellCpv
(
    const scalar T,
    const label celli
) const
{
    if (this->alpha2()[celli] < residualAlpha_.value())
    {
        return thermo1_->cellCpv(T, celli);
    }
    else if (this->alpha1()[celli] < residualAlpha_.value())
    {
        return thermo2_->cellCpv(T, celli);
    }
    scalar alphaRho1 = this->alpha1()[celli]*rho1_[celli];
    scalar alphaRho2 = this->alpha2()[celli]*rho2_[celli];
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
    if (this->alpha2().boundaryField()[patchi][facei] < residualAlpha_.value())
    {
        return thermo1_->patchFaceCpv(T, patchi, facei);
    }
    if (this->alpha1().boundaryField()[patchi][facei] < residualAlpha_.value())
    {
        return thermo2_->patchFaceCpv(T, patchi, facei);
    }
    scalar alphaRho1 =
        this->alpha1().boundaryField()[patchi][facei]
       *rho1_.boundaryField()[patchi][facei];
    scalar alphaRho2 =
        this->alpha2().boundaryField()[patchi][facei]
       *rho2_.boundaryField()[patchi][facei];

    return
        (
            alphaRho1*thermo1_->patchFaceCpv(T, patchi, facei)
          + alphaRho2*thermo2_->patchFaceCpv(T, patchi, facei)
        )/(alphaRho1 + alphaRho2);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::W() const
{
    return rho_/(this->alpha1()*rho1_/thermo1_->W() + this->alpha2()*rho2_/thermo2_->W());
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::W
(
    const label patchi
) const
{
    return
        rho_.boundaryField()[patchi]
       /(
            this->alpha1().boundaryField()[patchi]
           *rho1_.boundaryField()[patchi]
           /thermo1_->W(patchi)
          + this->alpha2().boundaryField()[patchi]
           *rho2_.boundaryField()[patchi]
           /thermo2_->W(patchi)
        );
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellW(const label celli) const
{
    return
        rho_[celli]
       /(
            this->alpha1()[celli]*rho1_[celli]/thermo1_->cellW(celli)
          + this->alpha2()[celli]*rho2_[celli]/thermo2_->cellW(celli)
        );
}

// ************************************************************************* //

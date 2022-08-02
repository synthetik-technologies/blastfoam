/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2019-04-29 Jeff Heylmun:    Simplified model
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

#include "multiPhaseModel.H"
#include "phaseSystem.H"
#include "fvMatrix.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiPhaseModel, 0);
    addToRunTimeSelectionTable
    (
        phaseModel,
        multiPhaseModel,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiPhaseModel::multiPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    fluidPhaseModel
    (
        fluid,
        phaseName,
        index,
        3 // This is just used to force the multiphase thermo
    ),
    thermo_(dynamicCast<multiphaseFluidBlastThermo>(thermoPtr_())),
    alphas_(thermo_.volumeFractions()),
    rhos_(thermo_.rhos()),
    alphaRhos_(alphas_.size()),
    alphaPhis_(alphas_.size()),
    alphaRhoPhis_(alphas_.size())
{
    thermo_.setTotalVolumeFractionPtr(*this);

    //- Temporarily Store read density
    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        0.0
    );
    alphaRho_ = dimensionedScalar(dimDensity, 0.0);

    wordList phaseNames(alphas_.size());
    forAll(alphas_, phasei)
    {
        phaseNames[phasei] = alphas_[phasei].group();
        sumAlpha += alphas_[phasei];
        word phaseName = alphas_[phasei].group();
        alphaRhos_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaRho", phaseName),
                    fluid.mesh().time().timeName(),
                    fluid.mesh()
                ),
                alphas_[phasei]*rhos_[phasei],
                rhos_[phasei].boundaryField().types()
            )
        );
        alphaRho_ += alphaRhos_[phasei];
        alphaPhis_.set
        (
            phasei,
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaPhi", phaseName),
                    fluid.mesh().time().timeName(),
                    fluid.mesh()
                ),
                fluid.mesh(),
                dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0.0)
            )
        );
        alphaRhoPhis_.set
        (
            phasei,
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaRhoPhi", phaseName),
                    fluid.mesh().time().timeName(),
                    fluid.mesh()
                ),
                fluid.mesh(),
                dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0.0)
            )
        );
    }

    this->fluxScheme_->phases() = phaseNames;

    // Reset density to correct value
    volScalarField& alpha = *this;
    alpha = sumAlpha;
    alpha.correctBoundaryConditions();

    rho_ = alphaRho_/Foam::max(sumAlpha, residualAlpha());

    solveAlpha(true);

    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiPhaseModel::~multiPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiPhaseModel::solve()
{
    dimensionedScalar dT = rho_.time().deltaT();
    dynamicCast<volScalarField>(*this) = 0.0;
    forAll(alphas_, phasei)
    {
        volScalarField deltaAlpha
        (
            fvc::div(alphaPhis_[phasei]) - alphas_[phasei]*fvc::div(phi_)
        );
        this->storeAndBlendDelta(deltaAlpha);

        volScalarField deltaAlphaRho(fvc::div(alphaRhoPhis_[phasei]));
        this->storeAndBlendDelta(deltaAlphaRho);

        this->storeAndBlendOld(alphas_[phasei]);
        alphas_[phasei] -= dT*deltaAlpha;
        alphas_[phasei].correctBoundaryConditions();

        this->storeAndBlendOld(alphaRhos_[phasei]);
        alphaRhos_[phasei].storePrevIter();
        alphaRhos_[phasei] -= dT*deltaAlphaRho;
        alphaRhos_[phasei].correctBoundaryConditions();

        *this += alphas_[phasei];
    }

    thermoPtr_->solve();
    phaseModel::solve();
}


void Foam::multiPhaseModel::postUpdate()
{
    // Solve phase mass
    bool needUpdate = false;
    alphaRho_.storePrevIter();
    forAll(rhos_, phasei)
    {
        bool alphaRhoUpdate = false;
        volScalarField& alpha(alphas_[phasei]);
        if (needSolve(alpha.name()))
        {
            //- Solve momentum equation (implicit stresses)
            fvScalarMatrix alphaEqn
            (
                fvm::ddt(alpha) - fvc::ddt(alpha)
             ==
                models().source(alpha)
            );
            constraints().constrain(alphaEqn);
            alphaEqn.solve();
            constraints().constrain(alpha);

            alphaRhoUpdate = true;
        }

        volScalarField& rho(rhos_[phasei]);
        if (needSolve(rho.name()))
        {
            dimensionedScalar rAlpha
            (
                thermo_.thermo(phasei).residualAlpha()
            );
            //- Solve momentum equation (implicit stresses)
            fvScalarMatrix rhoEqn
            (
                fvm::ddt(alpha, rho) - fvc::ddt(alphaRhos_[phasei])
              + fvm::ddt(rAlpha, rho)
              - fvc::ddt(rAlpha, rho)
             ==
                models().source(alpha, rho)
            );
            constraints().constrain(rhoEqn);
            rhoEqn.solve();
            constraints().constrain(rho);

            alphaRhoUpdate = true;
        }
        if (alphaRhoUpdate)
        {
            alphaRhos_[phasei] = alpha*rho;
            needUpdate = true;
        }
    }
    if (needUpdate)
    {
        alphaRho_ = alphaRhos_[0];
        for (label phasei = 1; phasei < alphaRhos_.size(); phasei++)
        {
            alphaRho_ += alphaRhos_[phasei];
        }
        rho_ = alphaRho_/Foam::max(*this, residualAlpha());
    }
    phaseModel::postUpdate();
}


void Foam::multiPhaseModel::update()
{
    fluxScheme_->update
    (
        alphas_,
        rhos_,
        U_,
        e_,
        p_,
        speedOfSound(),
        phi_,
        alphaPhiPtr_(),
        alphaRhoPhi_,
        alphaPhis_,
        alphaRhoPhis_,
        alphaRhoUPhi_,
        alphaRhoEPhi_
    );

    phaseModel::update();
    thermoPtr_->update();
}


void Foam::multiPhaseModel::correctVolumeFraction()
{
    // find largest volume fraction and set to 1-sum
    scalarList rAlphas(alphas_.size());
    SortableList<scalar> alphas(alphas_.size());
    forAll(rAlphas, phasei)
    {
        rAlphas[phasei] = thermo_.thermo(phasei).residualAlpha().value();
    }
    forAll(*this, celli)
    {
        forAll(alphas_, phasei)
        {
            alphas_[phasei][celli] =
                Foam::max
                (
                    Foam::min
                    (
                        alphas_[phasei][celli],
                        (*this)[celli]
                    ),
                    0.0
                );
            alphas[phasei] = alphas_[phasei][celli];
        }
        alphas.reverseSort();

        const label fixedPhase = alphas.indices()[0];

        scalar sumAlpha = 0.0;
        for (label phasei = 1; phasei < alphas.size(); phasei++)
        {
            // Conserve mass
            rhos_[phasei][celli] =
                Foam::max(alphaRhos_[phasei][celli], 0.0)
               /Foam::max(alphas[phasei], rAlphas[phasei]);

            sumAlpha += alphas[phasei];
        }
        sumAlpha = Foam::min(sumAlpha, (*this)[celli]);
        alphas_[fixedPhase][celli] = (*this)[celli] - sumAlpha;

        // Conserve mass
        rhos_[fixedPhase][celli] =
            Foam::max(alphaRhos_[fixedPhase][celli], 0.0)
           /Foam::max(alphas_[fixedPhase][celli], rAlphas[fixedPhase]);
    }

    forAll(alphas_, phasei)
    {
        alphas_[phasei].correctBoundaryConditions();
    }
}


void Foam::multiPhaseModel::decode()
{
    // Calculate densities
    alphaRho_ = dimensionedScalar("0", dimDensity, 0.0);

    forAll(alphas_, phasei)
    {
        alphaRhos_[phasei].max(0);
        rhos_[phasei] =
            alphaRhos_[phasei]
           /Foam::max
            (
                alphas_[phasei],
                thermo_.thermo(phasei).residualAlpha()
            );
        rhos_[phasei].correctBoundaryConditions();

        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];

        alphaRho_ += alphaRhos_[phasei];
    }
    volScalarField& alpha(*this);
    this->correctBoundaryConditions();

    rho_ = alphaRho_/Foam::max(alpha, residualAlpha());

    volScalarField alphaRhoLimited(alphaRho_);
    alphaRhoLimited.max(1e-10);
    U_.ref() = alphaRhoU_()/(alphaRhoLimited());
    U_.correctBoundaryConditions();

    alphaRhoU_.boundaryFieldRef() ==
        alphaRho_.boundaryField()*U_.boundaryField();

    e_.ref() = alphaRhoE_()/alphaRhoLimited() - 0.5*magSqr(U_());
    e_.correctBoundaryConditions();

    thermoPtr_->correct();
    thermoPtr_->speedOfSound() *= pos(alpha - residualAlpha());
    thermoPtr_->speedOfSound().max(small);

    // Update total energy because e may have changed
    alphaRhoE_ == alphaRho_*(e_ + 0.5*magSqr(U_));
}


void Foam::multiPhaseModel::encode()
{
    //- Scale volume fractions to new value
    alphaRho_ = dimensionedScalar(dimDensity, 0.0);
    volScalarField& alpha(*this);
    alpha = 0.0;
    forAll(alphas_, phasei)
    {
        alpha += alphas_[phasei];
        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
        alphaRho_ += alphaRhos_[phasei];
    }

    alphaRhoU_ == alphaRho_*U_;
    alphaRhoE_ == alphaRho_*(e_ + 0.5*magSqr(U_));
}


// ************************************************************************* //

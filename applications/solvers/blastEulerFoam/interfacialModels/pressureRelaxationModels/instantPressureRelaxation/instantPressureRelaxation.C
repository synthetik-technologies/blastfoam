/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.
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

#include "instantPressureRelaxation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(instantPressureRelaxation, 0);
    addToRunTimeSelectionTable
    (
        pressureRelaxationSolver,
        instantPressureRelaxation,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::instantPressureRelaxation::instantPressureRelaxation
(
    phaseSystem& fluid,
    interfacialPressureModelTable& interfacialPressureModels,
    pressureRelaxationModelTable& pressureRelaxationModels
)
:
    pressureRelaxationSolver
    (
        fluid,
        interfacialPressureModels,
        pressureRelaxationModels,
        false,
        false
    ),
    MultivariateEquation<scalar>
    (
        phaseModels_.size() + 1,
        scalarField(phaseModels_.size() + 1, 0.0),
        scalarField(phaseModels_.size() + 1, great)
    ),
    dict_(fluid.subDict("pressureSolverCoeffs")),
    rootSolver_(rootSolver::New(*this, dict_))
{
    pressureRelaxationSolver::nEqns_ = phaseModels_.size() + 1;
    rho0_.resize(phaseModels_.size());
    e0_.resize(phaseModels_.size());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::instantPressureRelaxation::~instantPressureRelaxation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::instantPressureRelaxation::f
(
     const scalarField& rhoPI,
    const label li,
    scalarField& fx
) const
{
    fx = scalarField(nEqns(), 0.0);
    scalar sumAlpha = 0;

    forAll(phaseModels_, phasei)
    {
        scalar alphaRho = phaseModels_[phasei].alphaRho()[li];
        if (alphaRho < 1e-10)
        {
            continue;
        }
        scalar alpha = alphaRho/max(rhoPI[phasei], 1e-10);
        thermos_[phasei].rho()[li] = rhoPI[phasei];
        scalar e = thermos_[phasei].calcCelle(rhoPI.last(), li);
        thermos_[phasei].he()[li] = e;

        fx[phasei] =
            2.0*rhoPI[phasei]*rho0_[phasei]*(e - e0_[phasei])
          + (rhoPI.last() + PI0_)*(rhoPI[phasei] - rho0_[phasei]);
        sumAlpha += alpha;
    }
    fx.last() = sumAlpha - 1.0;
}


void Foam::instantPressureRelaxation::jacobian
(
    const scalarField& rhoPI,
    const label li,
    scalarField& fx,
    RectangularMatrix<scalar>& J
) const
{
    J = scalarSquareMatrix(nEqns(), 0.0);
    fx = scalarField(nEqns(), 0.0);
    scalar sumAlpha = 0;
    Info<<rhoPI<<endl;
    forAll(phaseModels_, phasei)
    {
        scalar alphaRho = phaseModels_[phasei].alphaRho()[li];
        if (alphaRho < 1e-10)
        {
            continue;
        }
        scalar alpha = alphaRho/max(rhoPI[phasei], 1e-10);
        thermos_[phasei].rho()[li] = rhoPI[phasei];
        scalar e = thermos_[phasei].calcCelle(rhoPI.last(), li);
        thermos_[phasei].he()[li] = e;

        fx[phasei] =
            2.0*rhoPI[phasei]*rho0_[phasei]*(e - e0_[phasei])
          + (rhoPI.last() + PI0_)*(rhoPI[phasei] - rho0_[phasei]);

        J(phasei, phasei) =
            2.0*rho0_[phasei]
           *(
                (e - e0_[phasei])
              + rhoPI[phasei]*thermos_[phasei].celldpde(li)
            );
        J(phasei, phaseModels_.size()) = -alphaRho/sqr(rhoPI[phasei]);
        J(phaseModels_.size(), phasei) = rho0_[phasei] - rhoPI[phasei];

        sumAlpha += alpha;
    }
    fx.last() = sumAlpha - 1.0;
    J(phaseModels_.size(), phaseModels_.size()) = small;
}


Foam::scalar Foam::instantPressureRelaxation::solve
(
    const scalar& deltaT
)
{
    scalarField rhoPI(phaseModels_.size() + 1);
    forAll(phaseModels_[0], celli)
    {
        forAll(phaseModels_, phasei)
        {
            rho0_[phasei] = phaseModels_[phasei].rho()[celli];
            e0_[phasei] = phaseModels_[phasei].he()[celli];
            rhoPI[phasei] = rho0_[phasei];
        }
        PI0_ = fluid_.PI()[celli];
        rhoPI.last() = PI0_;

        rhoPI = rootSolver_->solve(rhoPI, celli);
    }
    forAll(thermos_, phasei)
    {
        thermos_[phasei].correct();
    }
    return deltaT;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
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

#include "plasticSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(plasticSolver, 0);
    const scalar plasticSolver::sqrt2By3 = sqrt(2.0/3.0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plasticSolver::plasticSolver
(
    const dictionary& dict
)
:
    maxDeltaErr_(dict.lookupOrDefault<scalar>("maxDeltaError", 0.01)),
    tolerance_(dict.lookupOrDefault<scalar>("plasticTolerance", 1e-8)),
    maxIter_(dict.lookupOrDefault<label>("maxIter", 100)),
    finiteDifference_(dict.lookupOrDefault<scalar>("finiteDifference", 1e-6))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::plasticSolver::~plasticSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::plasticSolver::fY
(
    const yieldStressModel& ys,
    const scalar epsilonPEq0,
    const scalar magSTrial,
    const scalar DLambda,
    const scalar muBar,
    const scalar J
) const
{
    return
        magSTrial
      - 2.0*muBar*DLambda
      - sqrt2By3
       *ys.currentSigmaY
        (
            epsilonPEq0,
            epsilonPEq0 + sqrt2By3*DLambda,
            J
        );
}

void Foam::plasticSolver::newtonLoop
(
    const yieldStressModel& ys,
    scalar& DLambda,
    scalar& sigmaY,
    const scalar epsilonPEq0,
    const scalar magSTrial,
    const scalar muBar,
    const scalar maxMagDEpsilon,
    const scalar J
) const
{
    label iter = 0;
    scalar residual = 1.0;

    do
    {
        const scalar fTrial =
            fY
            (
                ys,
                epsilonPEq0,
                magSTrial,
                DLambda,
                muBar,
                J
            );
        const scalar fTrialStep =
            fY
            (
                ys,
                epsilonPEq0,
                magSTrial,
                DLambda + finiteDifference_,
                muBar,
                J
            );
        const scalar dfTrialdLambda =
            (fTrialStep - fTrial)/finiteDifference_;

        if (mag(dfTrialdLambda) < small)
        {
            break;
        }
        residual = fTrial/dfTrialdLambda;
        DLambda -= residual;

        residual /= maxMagDEpsilon;
    } while (iter++ < maxIter_ && mag(residual) < tolerance_);

    if (iter == maxIter_)
    {
        WarningInFunction
            << "Plasticity Newton loop did not converge in "
            << maxIter_ << " iterations" << endl;
    }

    sigmaY = ys.currentSigmaY
    (
        epsilonPEq0,
        epsilonPEq0 + sqrt2By3*DLambda,
        J
    );
}

void Foam::plasticSolver::updatePlasticity
(
    const yieldStressModel& ys,
    symmTensor& plasticN,
    scalar& DLambda,
    scalar& sigmaY,
    const scalar sigmaY0,
    const scalar fTrial,
    const symmTensor& sTrial,
    const scalar epsilonPEq0,
    const scalar muBar,
    const scalar maxMagDEpsilon,
    const scalar J
) const
{
    if (fTrial < small)
    {
        plasticN = symmTensor::I;
        DLambda = 0.0;
        return;
    }

    const scalar magS = mag(sTrial);
    if (magS > small)
    {
        plasticN = sTrial/magS;
    }
    else
    {
        plasticN = symmTensor::I;
    }

    newtonLoop
    (
        ys,
        DLambda,
        sigmaY,
        epsilonPEq0,
        magS,
        muBar,
        maxMagDEpsilon,
        J
    );
}


// ************************************************************************* //


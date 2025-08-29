/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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

#include "NewtonRaphsonRootSolver.H"
#include "SVD.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace rootSolvers
{
    defineTypeNameAndDebug(NewtonRaphson, 0);
    addToRunTimeSelectionTable
    (
        rootSolver,
        NewtonRaphson,
        dictionaryOne
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rootSolvers::NewtonRaphson::NewtonRaphson
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
:
    rootSolver(eqn, dict),
    beta_(dict.lookupOrDefault<scalar>("beta", 1.0))
{}


Foam::rootSolvers::NewtonRaphson::NewtonRaphson
(
    const scalarMultivariateEquation& eqn,
    const NewtonRaphson& solver
)
:
    rootSolver(eqn, solver),
    beta_(solver.beta_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rootSolvers::NewtonRaphson::~NewtonRaphson()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::rootSolvers::NewtonRaphson::findRoots
(
    const scalarList& x0,
    const scalarList& xLow,
    const scalarList& xHigh,
    const label li
) const
{
    initialise(x0);
    scalarField xOld(x0);
    tmp<scalarField> xNewTmp(new scalarField(x0));
    scalarField& xNew = xNewTmp.ref();
    scalarField f(xNew.size());
    RectangularMatrix<scalar> J(xNew.size());
    eqns_.jacobian(xOld, li, f, J);

    const scalarList& lower = eqns_.lowerLimits();
    const scalarList& upper = eqns_.upperLimits();

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        scalarField delta(-(SVDinv(J)*f));

        // Relax delta
        if (beta_ != 1.0)
        {
            delta *= 1.0/(1.0 + beta_*sum(magSqr(delta)));
        }

        xNew = xOld + delta;

        scalar fac = 0.5;
        // eqns_.limit(xNew);
        forAll(xNew, i)
        {
            if (xNew[i] < lower[i])
            {
                xNew[i] = (1.0 - fac)*xOld[i] + fac*lower[i];
            }
            else if (xNew[i] > upper[i])
            {
                xNew[i] = (1.0 - fac)*xOld[i] + fac*upper[i];
            }
        }

        if (converged(xOld, xNew))
        {
            break;
        }

        eqns_.jacobian(xNew, li, f, J);

        printStepInformation(xNew);

        xOld = xNew;
    }
    printFinalInformation(xNew);

    return xNewTmp;
}

// ************************************************************************* //

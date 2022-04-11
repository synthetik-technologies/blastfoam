/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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

#include "goodBroydenRootSolver.H"
#include "SVD.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(goodBroydenRootSolver, 0);
    addToRunTimeSelectionTable
    (
        rootSolver,
        goodBroydenRootSolver,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        rootSolver,
        goodBroydenRootSolver,
        dictionaryOne
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::goodBroydenRootSolver::goodBroydenRootSolver
(
    const scalarMultivariateEquation& eqns,
    const dictionary& dict
)
:
    rootSolver(eqns, dict)
{
    if (dict.found("dx"))
    {
        eqns_.setDX(dict.lookup<scalarField>("dx"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::goodBroydenRootSolver::~goodBroydenRootSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::goodBroydenRootSolver::findRoots
(
    const scalarList& x0,
    const scalarList& xLow,
    const scalarList& xHigh,
    const label li
) const
{
    initialise(x0);
    tmp<scalarField> xTmp(new scalarField(x0));
    scalarField& x = xTmp.ref();

    scalarField xOld(x0);
    scalarRectangularMatrix J(eqns_.nEqns(), x0.size());
    scalarField fOld(eqns_.nEqns());
    eqns_.jacobian(x0, li, fOld, J);
    scalarField f(fOld);
    scalarRectangularMatrix Jinv(SVDinv(J));

    scalarRectangularMatrix dx(x.size(), 1);
    scalarRectangularMatrix df(x.size(), 1);

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        // Save previous time step
        xOld = x;
        fOld = f;

        // Increment by the jacobian
        x -= Jinv*f;

        // Limit to the bounds to the equation
        eqns_.limit(x);

        // update f
        eqns_.FX(x, li, f);

        // Check for convergence
        if (converged(xOld, x, f))
        {
            break;
        }

        // Update changes in x and f
        forAll(dx, i)
        {
            dx(i, 0) = x[i] - xOld[i];
            df(i, 0) = f[i] - fOld[i];
        }

        // Update Jinv
        Jinv =
            Jinv
          + ((dx - Jinv*df)*(dx.T()*Jinv))
           /stabilise(((dx.T()*Jinv)*df)(0, 0), small);

        printStepInformation(x);
    }
    printFinalInformation(x);

    return xTmp;
}

// ************************************************************************* //

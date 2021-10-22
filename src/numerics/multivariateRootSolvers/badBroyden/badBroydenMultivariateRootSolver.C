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

#include "badBroydenMultivariateRootSolver.H"
#include "SVD.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(badBroydenMultivariateRootSolver, 0);
    addToRunTimeSelectionTable
    (
        multivariateRootSolver,
        badBroydenMultivariateRootSolver,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        multivariateRootSolver,
        badBroydenMultivariateRootSolver,
        dictionaryOne
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::badBroydenMultivariateRootSolver::badBroydenMultivariateRootSolver
(
    const scalarMultivariateEquation& eqns,
    const dictionary& dict
)
:
    multivariateRootSolver(eqns, dict),
    dx_(dict.lookupOrDefault("dx", scalarField(eqns.nEqns(), 1e-3)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::badBroydenMultivariateRootSolver::findRoots
(
    const scalarList& x0,
    const scalarList& xLow,
    const scalarList& xHigh,
    const label li
) const
{
    scalarField xOld(x0);
    scalarRectangularMatrix J(eqns_.calculateJacobian(x0, dx_, li));

    tmp<scalarField> xTmp(new scalarField(x0));
    scalarField& x = xTmp.ref();
    scalarRectangularMatrix dx(x.size(), 1);
    scalarRectangularMatrix df(x.size(), 1);

    scalarRectangularMatrix Jinv(SVDinv(J));

    scalarField fOld(eqns_.f(x0, li));
    scalarField f(fOld);

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        xOld = x;
        fOld = f;

        scalarField delta(-Jinv*f);
        x += delta;

        eqns_.limit(x);

        if (converged(delta))
        {
            break;
        }

        // update f
        eqns_.f(x, li, f);

        forAll(dx, i)
        {
            dx(i, 0) = x[i] - xOld[i];
            df(i, 0) = f[i] - fOld[i];
        }

        // Update Jinv
        Jinv =
            Jinv
          + ((dx - Jinv*df)*df.T())/stabilise((df.T()*df)(0, 0), small);

        printStepInformation(x);
    }
    printFinalInformation();

    return xTmp;
}

// ************************************************************************* //

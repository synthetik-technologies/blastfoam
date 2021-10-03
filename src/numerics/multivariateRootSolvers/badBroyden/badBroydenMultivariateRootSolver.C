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

Foam::tmp<Foam::scalarField> Foam::badBroydenMultivariateRootSolver::solve
(
    const scalarField& x0,
    const label li
) const
{
    scalarField xOld(x0);
    scalarSquareMatrix J(x0.size());
    forAll(x0, cmptj)
    {
        scalarField x1(x0);
        scalarField x2(x0);
        x1[cmptj] -= dx_[cmptj];
        x2[cmptj] += dx_[cmptj];

        scalarField f1(eqns_.f(x1, li));
        scalarField f2(eqns_.f(x2, li));

        forAll(x0, cmpti)
        {
            J(cmpti, cmptj) =
                (f2[cmpti] - f1[cmpti])/(x2[cmptj] - x1[cmptj]);
        }
    }

    tmp<scalarField> xTmp(new scalarField(x0));
    scalarField& x = xTmp.ref();
    scalarRectangularMatrix dx(x.size(), 1);
    scalarRectangularMatrix df(x.size(), 1);

    scalarSquareMatrix Jinv(SVDinv(J));

    scalarField fOld(eqns_.f(x0, li));
    scalarField f(fOld);

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        xOld = x;
        fOld = f;

        scalarField delta(-Jinv*f);
        x += delta;

        if (converged(delta))
        {
            return xTmp;
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
          + scalarSquareMatrix
            (
                ((dx - Jinv*df)*df.T())/(df.T()*df)(0, 0)
            );
    }
    WarningInFunction
        << "Could not converge to the given multivariateMultivariateRoot." << endl;

    return xTmp;
}

// ************************************************************************* //

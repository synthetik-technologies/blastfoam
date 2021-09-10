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

#include "HalleyRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HalleyRootSolver, 0);
    addToRunTimeSelectionTable(rootSolver, HalleyRootSolver, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HalleyRootSolver::HalleyRootSolver
(
    const rootSystem& eqn,
    const dictionary& dict
)
:
    rootSolver(eqn, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::HalleyRootSolver::solve
(
    const scalar x0,
    const label li
) const
{
    scalar xOld = x0;
    scalar xNew = x0;
    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        scalar f = eqn_.f(xOld, li);
        scalar fp = eqn_.dfdx(xOld, li);
        scalar fpp = eqn_.d2fdx2(xOld, li);

        xNew = xOld - 2.0*f*fp/stabilise(2.0*sqr(fp) - f*fpp, tolerance_);

        if (converged(xNew - xOld))
        {
            return xNew;
        }
        xOld = xNew;

    }
    WarningInFunction
        << "Could not converge to the given root." << endl;

    return xNew;
}

// ************************************************************************* //

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

#include "NewtonRaphsonMultivariateRootSolver.H"
#include "SVD.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(NewtonRaphsonMultivariateRootSolver, 0);
    addToRunTimeSelectionTable
    (
        multivariateRootSolver,
        NewtonRaphsonMultivariateRootSolver,
        dictionaryOne
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NewtonRaphsonMultivariateRootSolver::NewtonRaphsonMultivariateRootSolver
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
:
    multivariateRootSolver(eqn, dict),
    beta_(dict.lookupOrDefault<scalar>("beta", 1.0))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::NewtonRaphsonMultivariateRootSolver::solve
(
    const scalarField& x0,
    const label li
) const
{
    scalarField xOld(x0);
    tmp<scalarField> xNewTmp(new scalarField(x0));
    scalarField& xNew = xNewTmp.ref();
    scalarField f(xNew.size());
    scalarSquareMatrix J(xNew.size());

    eqns_.jacobian(xOld, li, f, J);

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        scalarField delta(-(SVDinv(J)*f));

        if (converged(delta))
        {
            return xNewTmp;
        }
        // Relax delta
        if (beta_ != 1.0)
        {
            delta *= 1.0/(1.0 + beta_*sum(magSqr(delta)));
        }

        xNew = xOld + delta;
        limit(xNew);
        xOld = xNew;
        eqns_.jacobian(xOld, li, f, J);
    }
    WarningInFunction
        << "Could not converge. Final error=" << error_ << endl;

    return xNewTmp;
}

// ************************************************************************* //

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

#include "univariateRootSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(univariateRootSolver, 0);
    defineRunTimeSelectionTable(univariateRootSolver, dictionaryZero);
    defineRunTimeSelectionTable(univariateRootSolver, dictionaryOne);
    defineRunTimeSelectionTable(univariateRootSolver, dictionaryTwo);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::univariateRootSolver::converged(const scalar error) const
{
    errors_[0] = mag(error);
    return errors_[0] < tolerances_[0];
}



void Foam::univariateRootSolver::printStepInformation(const scalar val) const
{
    if (debug)
    {
        Info<< "Step " << stepi_
            << ", error=" << errors_[0]
            << ", value: " << val << endl;
    }
}


Foam::scalar
Foam::univariateRootSolver::printFinalInformation(const scalar val) const
{
    if (stepi_ < maxSteps_ && debug > 1)
    {
        Info<< "Converged in " << stepi_ << " iterations"
            << ", final error=" << errors_[0];
    }
    else if (stepi_ >= maxSteps_ && debug)
    {
        WarningInFunction
            << "Did not converge, final error= " << errors_[0] << endl;
    }
    return val;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateRootSolver::univariateRootSolver
(
    const multivariateEquation<scalar>& eqn,
    const dictionary& dict
)
:
    rootSolver(eqn, dict),
    eqn_(dynamicCast<const equation>(eqn))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::univariateRootSolver::~univariateRootSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::univariateRootSolver::solve(const scalar x0) const
{
    return solve(x0, eqn_.lower(), eqn_.upper(), 0);
}


Foam::scalar Foam::univariateRootSolver::solve
(
    const scalar x0,
    const label li
) const
{
    return solve(x0, eqn_.lower(), eqn_.upper(), li);
}


Foam::scalar Foam::univariateRootSolver::solve
(
    const scalar x0,
    const scalar xLow,
    const scalar xHigh
) const
{
    return solve(x0, xLow, xHigh, 0);
}


Foam::scalar Foam::univariateRootSolver::solve
(
    const scalar x0,
    const scalar xLow,
    const scalar xHigh,
    const label li
) const
{
    return this->findRoot(x0, xLow, xHigh, li);
}


Foam::tmp<Foam::scalarField> Foam::univariateRootSolver::findRoots
(
    const scalarField& x0,
    const scalarField& xLow,
    const scalarField& xHigh,
    const label li
) const
{
    return tmp<scalarField>
    (
        new scalarField(1, this->findRoot(x0[0], xLow[0], xHigh[0], li))
    );
}


// ************************************************************************* //

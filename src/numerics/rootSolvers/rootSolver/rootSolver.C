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

#include "rootSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rootSolver, 0);
    defineRunTimeSelectionTable(rootSolver, dictionaryZero);
    defineRunTimeSelectionTable(rootSolver, dictionaryOne);
    defineRunTimeSelectionTable(rootSolver, dictionaryTwo);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::rootSolver::converged(const scalar error) const
{
    error_ = mag(error);
    return error_ < tolerance_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rootSolver::rootSolver(const scalarEquation& eqn, const dictionary& dict)
:
    eqn_(eqn),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", 1e-6)),
    maxSteps_(dict.lookupOrDefault<scalar>("maxSteps", 100)),
    stepi_(0),
    error_(great)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::rootSolver::solve() const
{
    return this->solve
    (
        (eqn_.lower() + eqn_.upper())*0.5,
        eqn_.lower(),
        eqn_.upper(),
        0
    );
}


Foam::scalar Foam::rootSolver::solve(const scalar x0) const
{
    return this->solve(x0, eqn_.lower(), eqn_.upper(), 0);
}


Foam::scalar Foam::rootSolver::solve(const scalar x0, const label li) const
{
    return this->solve(x0, eqn_.lower(), eqn_.upper(), li);
}


Foam::scalar Foam::rootSolver::solve
(
    const scalar x0,
    const scalar xLow,
    const scalar xHigh
) const
{
    return this->solve(x0, xLow, xHigh, 0);
}


// ************************************************************************* //

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::univariateRootSolver> Foam::univariateRootSolver::New
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
{
    word rootSolverTypeName(dict.lookup("solver"));
    label nDeriv = eqn.nDerivatives();
    Info<< "Selecting root solver " << rootSolverTypeName << endl;
    if (debug)
    {
        Info<< "    detected " << nDeriv << " implemented derivatives" << endl;
    }

    if (nDeriv <= 0)
    {
        dictionaryZeroConstructorTable::iterator cstrIter =
            dictionaryZeroConstructorTablePtr_->find(rootSolverTypeName);

        if (cstrIter == dictionaryZeroConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown univariateRootSolver type "
                << rootSolverTypeName << nl << nl
                << "Valid univariateRootSolver for no derivatives are : " << endl
                << dictionaryZeroConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }
        return autoPtr<univariateRootSolver>(cstrIter()(eqn, dict));
    }
    else if (nDeriv == 1)
    {
        dictionaryOneConstructorTable::iterator cstrIter =
            dictionaryOneConstructorTablePtr_->find(rootSolverTypeName);

        if (cstrIter == dictionaryOneConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown univariateRootSolver type "
                << rootSolverTypeName << nl << nl
                << "Valid univariateRootSolver for one derivative are : " << endl
                << dictionaryOneConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }
        return autoPtr<univariateRootSolver>(cstrIter()(eqn, dict));
    }
    dictionaryTwoConstructorTable::iterator cstrIter =
        dictionaryTwoConstructorTablePtr_->find(rootSolverTypeName);

    if (cstrIter == dictionaryTwoConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown univariateRootSolver type "
            << rootSolverTypeName << nl << nl
            << "Valid univariateRootSolver for are : " << endl
            << dictionaryTwoConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return autoPtr<univariateRootSolver>(cstrIter()(eqn, dict));
}


Foam::autoPtr<Foam::rootSolver> Foam::univariateRootSolver::NewUnivariate
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
{
    autoPtr<univariateRootSolver> uRootSolver(New(eqn, dict));
    return autoPtr<rootSolver>(uRootSolver.ptr());
}


// ************************************************************************* //

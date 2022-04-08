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
    const scalarEquation& eqn,
    const dictionary& dict
)
{
    return New(dict.lookup<word>("solver"), eqn, dict);
}


Foam::autoPtr<Foam::univariateRootSolver> Foam::univariateRootSolver::New
(
    const word& rootSolverType,
    const scalarEquation& eqn,
    const dictionary& dict
)
{
    label nDeriv = eqn.nDerivatives();
    Info<< "Selecting univariate root solver: " << rootSolverType << endl;
    if (nDeriv <= 0)
    {
        dictionaryZeroConstructorTable::iterator cstrIter =
            dictionaryZeroConstructorTablePtr_->find(rootSolverType);

        if (cstrIter == dictionaryZeroConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown univariateRootSolver type "
                << rootSolverType << nl << nl
                << "Valid univariateRootSolver for no derivatives are : " << endl
                << dictionaryZeroConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }
        return autoPtr<univariateRootSolver>(cstrIter()(eqn, dict));
    }
    else if (nDeriv == 1)
    {
        dictionaryOneConstructorTable::iterator cstrIter =
            dictionaryOneConstructorTablePtr_->find(rootSolverType);

        if (cstrIter == dictionaryOneConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown univariateRootSolver type "
                << rootSolverType << nl << nl
                << "Valid univariateRootSolver for one derivative are : " << endl
                << dictionaryOneConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }
        return autoPtr<univariateRootSolver>(cstrIter()(eqn, dict));
    }
    else
    {
        dictionaryTwoConstructorTable::iterator cstrIter =
            dictionaryTwoConstructorTablePtr_->find(rootSolverType);

        if (cstrIter == dictionaryTwoConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown univariateRootSolver type "
                << rootSolverType << nl << nl
                << "Valid univariateRootSolver for are : " << endl
                << dictionaryTwoConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }
        return autoPtr<univariateRootSolver>(cstrIter()(eqn, dict));
    }
}


// ************************************************************************* //

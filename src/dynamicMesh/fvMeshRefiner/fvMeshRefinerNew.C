/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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


#include "fvMeshHexRefiner.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fvMeshRefiner> Foam::fvMeshRefiner::New
(
    fvMesh& mesh
)
{
    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh.time().constant(),
            mesh.dbDir(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word refinerType = fvMeshHexRefiner::typeName;
    if (dynamicMeshDict.found("refiner") || mesh.nSolutionD() < 2)
    {
        refinerType = dynamicMeshDict.lookup<word>("refiner");
    }

    Info<< "Selecting fvMeshRefiner " << refinerType << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(refinerType);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fvMeshRefiner type "
            << refinerType << nl << nl
            << "Valid fvMeshRefiner are :" << endl
            << fvMeshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<fvMeshRefiner>(cstrIter()(mesh));
}


Foam::autoPtr<Foam::fvMeshRefiner> Foam::fvMeshRefiner::New
(
    fvMesh& mesh,
    const dictionary& dict,
    const bool force,
    const bool read
)
{
    word refinerType = fvMeshHexRefiner::typeName;
    if (dict.found("refiner") || mesh.nSolutionD() < 2)
    {
        refinerType = dict.lookup<word>("refiner");
    }

    Info<< "Selecting fvMeshRefiner " << refinerType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(refinerType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fvMeshRefiner type "
            << refinerType << nl << nl
            << "Valid fvMeshRefiner are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<fvMeshRefiner>(cstrIter()(mesh, dict, force, read));
}


// ************************************************************************* //

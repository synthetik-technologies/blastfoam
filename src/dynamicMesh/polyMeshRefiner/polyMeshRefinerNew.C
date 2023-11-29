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


#include "polyMeshRefiner.H"
#include "polyMeshHexRefiner.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyMeshRefiner> Foam::polyMeshRefiner::New
(
    polyMesh& mesh
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

    word refinerType = polyMeshHexRefiner::typeName;
    if (dynamicMeshDict.found("refiner") || mesh.nSolutionD() < 2)
    {
        refinerType = dynamicMeshDict.lookup<word>("refiner");
    }

    return New(refinerType, mesh);
}

Foam::autoPtr<Foam::polyMeshRefiner> Foam::polyMeshRefiner::New
(
    const word& refinerType,
    polyMesh& mesh
)
{
    polyMeshConstructorTable::iterator cstrIter =
        polyMeshConstructorTablePtr_->find(refinerType);

    if (cstrIter == polyMeshConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown polyMeshRefiner type "
            << refinerType << nl << nl
            << "Valid polyMeshRefiner are :" << endl
            << polyMeshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<polyMeshRefiner>(cstrIter()(mesh));
}


Foam::autoPtr<Foam::polyMeshRefiner> Foam::polyMeshRefiner::New
(
    polyMesh& mesh,
    const dictionary& dict,
    const bool force,
    const bool read
)
{
    word refinerType = polyMeshHexRefiner::typeName;
    if (dict.found("refiner") || mesh.nSolutionD() < 2)
    {
        refinerType = dict.lookup<word>("refiner");
    }

    return New(refinerType, mesh, dict, force, read);
}


Foam::autoPtr<Foam::polyMeshRefiner> Foam::polyMeshRefiner::New
(
    const word& refinerType,
    polyMesh& mesh,
    const dictionary& dict,
    const bool force,
    const bool read
)
{
    Info<< "Selecting polyMeshRefiner " << refinerType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(refinerType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown polyMeshRefiner type "
            << refinerType << nl << nl
            << "Valid polyMeshRefiner are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<polyMeshRefiner>(cstrIter()(mesh, dict, force, read));
}

// ************************************************************************* //

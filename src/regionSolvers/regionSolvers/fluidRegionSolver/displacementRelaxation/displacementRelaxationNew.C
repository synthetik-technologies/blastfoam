/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
------------------------------------------------------------------------*/

#include "displacementRelaxation.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::displacementRelaxation> Foam::displacementRelaxation::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    word relaxType("fixed");
    if (dict.isDict(mesh.name())&& dict.subDict(mesh.name()).found("relaxation"))
    {
        dict.subDict(mesh.name()).lookup("relaxation") >> relaxType;
    }
    else if (dict.found("relaxation"))
    {
        dict.lookup("relaxation") >> relaxType;
    }

    Info<< "Selecting relaxation method: " << relaxType
        << " for region " << mesh.name() << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(relaxType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown relaxation type "
            << relaxType << endl << endl
            << "Valid regionSolvers are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mesh, dict);
}


// ************************************************************************* //

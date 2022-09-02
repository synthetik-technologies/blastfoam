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

\*---------------------------------------------------------------------------*/

#include "fluxScheme.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluxScheme> Foam::fluxScheme::NewSingle
(
    const fvMesh& mesh
)
{
    word fluxSchemeType(mesh.schemesDict().lookup<word>("fluxScheme"));

    Info<< "Selecting fluxScheme: " << fluxSchemeType << endl;

    singlePhaseConstructorTable::iterator cstrIter =
        singlePhaseConstructorTablePtr_->find(fluxSchemeType);

    if (cstrIter == singlePhaseConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown singlePhase fluxScheme type "
            << fluxSchemeType << endl << endl
            << "Valid fluxScheme types are : " << endl
            << singlePhaseConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mesh);
}


Foam::autoPtr<Foam::fluxScheme> Foam::fluxScheme::NewMulti
(
    const fvMesh& mesh
)
{
    word fluxSchemeType(mesh.schemesDict().lookup<word>("fluxScheme"));

    Info<< "Selecting fluxScheme: " << fluxSchemeType << endl;

    multiphaseConstructorTable::iterator cstrIter =
        multiphaseConstructorTablePtr_->find(fluxSchemeType);

    if (cstrIter == multiphaseConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown multiphase fluxScheme type "
            << fluxSchemeType << endl << endl
            << "Valid fluxScheme types are : " << endl
            << multiphaseConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mesh);
}


Foam::autoPtr<Foam::fluxScheme> Foam::fluxScheme::NewInterface
(
    const fvMesh& mesh
)
{
    word fluxSchemeType(mesh.schemesDict().lookup<word>("fluxScheme"));

    Info<< "Selecting fluxScheme: " << fluxSchemeType << endl;

    interfaceConstructorTable::iterator cstrIter =
        interfaceConstructorTablePtr_->find(fluxSchemeType);

    if (cstrIter == interfaceConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown interface fluxScheme type "
            << fluxSchemeType << endl << endl
            << "Valid fluxScheme types are : " << endl
            << interfaceConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mesh);
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "phaseFluxScheme.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseFluxScheme> Foam::phaseFluxScheme::New
(
    const fvMesh& mesh,
    const word& name
)
{
    word fluxSchemeType
    (
        mesh.schemesDict().subDict("fluxSchemes").subDict(name).lookup("fluxScheme")
    );

    Info<< "Selecting phaseFluxScheme: " << fluxSchemeType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(fluxSchemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown phaseFluxScheme type "
            << fluxSchemeType << endl << endl
            << "Valid phaseFluxScheme types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mesh, name);
}


Foam::autoPtr<Foam::phaseFluxScheme> Foam::phaseFluxScheme::NewSolid
(
    const fvMesh& mesh,
    const word& name
)
{
    word fluxSchemeType
    (
        mesh.schemesDict().subDict("fluxSchemes").subDict(name).lookup("fluxScheme")
    );

    Info<< "Selecting phaseFluxScheme: " << fluxSchemeType << endl;

    solidConstructorTable::iterator cstrIter =
        solidConstructorTablePtr_->find(fluxSchemeType);

    if (cstrIter == solidConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown phaseFluxScheme type "
            << fluxSchemeType << endl << endl
            << "Valid phaseFluxScheme types are : " << endl
            << solidConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mesh, name);
}


// ************************************************************************* //

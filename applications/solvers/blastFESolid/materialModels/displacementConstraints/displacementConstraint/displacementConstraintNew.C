/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
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

#include "displacementConstraint.H"

// * * * * * * * * * * * * * * Static Data Fucntions * * * * * * * * * * * * //

Foam::autoPtr<Foam::displacementConstraint>
Foam::displacementConstraint::New
(
    pointVectorField& D,
    pointVectorField& U,
    const dictionary& dict
)
{
    return New
    (
        dict.lookup<word>("type"),
        dict.dictName(),
        D,
        U,
        dict
    );
}


Foam::autoPtr<Foam::displacementConstraint>
Foam::displacementConstraint::New
(
    const word& type,
    const word& name,
    pointVectorField& D,
    pointVectorField& U,
    const dictionary& dict
)
{
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown displacement constraint " << type << nl
            << "Valid displacement constraint are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return cstrIter()(name, D, U, dict);
}


// ************************************************************************* //


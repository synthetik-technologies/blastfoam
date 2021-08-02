/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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

#include "solidBlastThermo.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidBlastThermo> Foam::solidBlastThermo::NewBasic
(
    const word& phaseName,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
{
     basicSolidConstructorTable::iterator cstrIter =
        lookupThermo<basicSolidConstructorTable>
        (
            dict,
            basicSolidConstructorTablePtr_
        );

    return cstrIter()(phaseName, mesh, dict, master, masterName);
}


Foam::autoPtr<Foam::solidBlastThermo> Foam::solidBlastThermo::NewDetonating
(
    const word& phaseName,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
{
    detonatingSolidConstructorTable::iterator cstrIter =
        lookupThermo<detonatingSolidConstructorTable>
        (
            dict,
            detonatingSolidConstructorTablePtr_
        );

    return cstrIter()(phaseName, mesh, dict, master, masterName);
}


Foam::autoPtr<Foam::solidBlastThermo> Foam::solidBlastThermo::New
(
    const word& phaseName,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
{
    const word type(dict.lookup<word>("type"));

    if (type == "basic")
    {
        return NewBasic(phaseName, mesh, dict, master, masterName);
    }
    else if (type == "detonating")
    {
        return NewDetonating(phaseName, mesh, dict, master, masterName);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown solidThermo type " << type << nl << nl
            << "Valid " << solidBlastThermo::typeName << " types are:" << nl
            << "basic" << nl
            << "detonating" << nl
            << abort(FatalError);
    }
    return autoPtr<solidBlastThermo>();
}

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2025-06-09 Jeff Heylmun     : Modified write function
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

#define list blastList
#include "fvMesh.H"

#include "blastListFvMeshTopoChanger.H"
#include "polyTopoChangeMap.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(blastList, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, blastList, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::blastList::blastList
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    fvMeshTopoChanger(mesh)
{
    const dictionary& solversDict = dict.subDict("topoChangers");

    forAllConstIter(dictionary, solversDict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& dict = iter().dict();

            list_.insert
            (
                name,
                fvMeshTopoChanger::New(mesh, dict).ptr()
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::blastList::~blastList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::blastList::update()
{
    bool updated = false;

    forAllIter(PtrDictionary<fvMeshTopoChanger>, list_, iter)
    {
        updated = iter().update() || updated;
        mesh().topoChanged_ = updated;
    }
    if (updated)
    {
        mesh().moving_ = false;
    }

    return updated;
}


void Foam::fvMeshTopoChangers::blastList::topoChange(const polyTopoChangeMap& map)
{
    forAllIter(PtrDictionary<fvMeshTopoChanger>, list_, iter)
    {
        iter().topoChange(map);
    }
}


void Foam::fvMeshTopoChangers::blastList::mapMesh(const polyMeshMap& map)
{
    forAllIter(PtrDictionary<fvMeshTopoChanger>, list_, iter)
    {
        iter().mapMesh(map);
    }
}


void Foam::fvMeshTopoChangers::blastList::distribute
(
    const polyDistributionMap& map
)
{
    forAllIter(PtrDictionary<fvMeshTopoChanger>, list_, iter)
    {
        iter().distribute(map);
    }
}


bool Foam::fvMeshTopoChangers::blastList::write(const bool w) const
{
    bool good = true;
    forAllConstIter(PtrDictionary<fvMeshTopoChanger>, list_, iter)
    {
        good = iter().write(w) && good;
    }
    return good;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License

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

#include "mappedPatchSelectorList.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedPatchSelectorList, 0);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedPatchSelectorList::mappedPatchSelectorList
(
    const polyMesh& mesh
)
:
    MappedPatchSelectorList(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedPatchSelectorList::~mappedPatchSelectorList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::mappedPatchSelectorList::movePoints()
{
    forAllIter
    (
        HashPtrTable<mappedPatchSelector>,
        patches_,
        iter
    )
    {
        iter()->clearOut();
    }
    return true;
}


void Foam::mappedPatchSelectorList::updateMesh(const mapPolyMesh& mpm)
{
    forAllIter
    (
        HashPtrTable<mappedPatchSelector>,
        patches_,
        iter
    )
    {
        iter()->clearOut();
    }
}



void Foam::mappedPatchSelectorList::reorderPatches
(
    const labelUList& newToOld,
    const bool validBoundary
)
{
    hashedWordList toc(patches_.toc());
    patches_.clear();
    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (mesh_.boundaryMesh().set(patchi))
        {
            const polyPatch& pp = mesh_.boundaryMesh()[patchi];
            if (toc.found(pp.name()))
            {
                patches_.insert(pp.name(), new mappedPatchSelector(pp));
            }
        }
    }
}


void Foam::mappedPatchSelectorList::addPatch(const label patchi)
{}


// * * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * //

const Foam::mappedPatchSelector&
Foam::mappedPatchSelectorList::operator()(const polyPatch& pp) const
{
    if (!patches_.found(pp.name()))
    {
        patches_.insert(pp.name(), new mappedPatchSelector(pp));
    }

    return *patches_[pp.name()];
}


const Foam::mappedPatchSelector&
Foam::mappedPatchSelectorList::operator()(const pointPatch& pp) const
{
    if (!patches_.found(pp.name()))
    {
        patches_.insert
        (
            pp.name(),
            new mappedPatchSelector
            (
                dynamicCast<const polyMesh>
                (
                    pp.boundaryMesh().mesh().thisDb()
                ).boundaryMesh()[pp.index()]
            )
        );
    }

    return *patches_[pp.name()];
}

// ************************************************************************* //

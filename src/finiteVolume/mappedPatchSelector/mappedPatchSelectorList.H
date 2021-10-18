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

Class
    Foam::mappedPatchSelectorList

Description
    Mesh object used to collect all mapped patches so that they can be cleared


Note

SourceFiles
    mappedPatchSelectorList.C

\*---------------------------------------------------------------------------*/

#ifndef mappedPatchSelectorList_H
#define mappedPatchSelectorList_H

#include "polyMesh.H"
#include "pointMesh.H"
#include "mappedPatchSelector.H"
#include "regIOobject.H"
#include "HashPtrTable.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mappedPatchSelectorList;
typedef MeshObject
<
    polyMesh,
    PatchMeshObject,
    mappedPatchSelectorList
> MappedPatchSelectorList;

/*---------------------------------------------------------------------------*\
                       Class mappedPatchSelectorList Declaration
\*---------------------------------------------------------------------------*/

class mappedPatchSelectorList
:
    public MappedPatchSelectorList
{
// Private member data

    // Table of mapped patch selectors
    mutable HashPtrTable<mappedPatchSelector> patches_;


public:

    //- Runtime type information
    TypeName("mappedPatchSelectorList");


    // Constructors

        //- Construct from patch
        mappedPatchSelectorList(const polyMesh&);


    //- Destructor
    virtual ~mappedPatchSelectorList();


    // Move points
    virtual bool movePoints();

    //- Update mesh
    virtual void updateMesh(const mapPolyMesh& mpm);

    //- Reorder patches
    virtual void reorderPatches
    (
        const labelUList& newToOld,
        const bool validBoundary
    );

    //- Inserted patch at patchi
    virtual void addPatch(const label patchi);

    // Operators

        //- Return the mappedPatchSelector for the given polyPatch
        const mappedPatchSelector& operator()(const polyPatch&) const;

        //- Return the mappedPatchSelector for the given pointPatch
        const mappedPatchSelector& operator()(const pointPatch&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

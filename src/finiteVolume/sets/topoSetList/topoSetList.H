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

Class
    Foam::topoSetList

Description
    Creation and storage of sets and zone

SourceFiles
    topoSetSource.C

\*---------------------------------------------------------------------------*/

#ifndef topoSetList_H
#define topoSetList_H

#include "fvMesh.H"
#include "HashPtrTable.H"
#include "topoSet.H"
#include "cellZoneSet.H"
#include "faceZoneSet.H"
#include "pointZoneSet.H"
#include "NamedEnum.H"
#include "RefineBalanceMeshObject.H"
#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class topoSetList;
typedef MeshObject
<
    fvMesh,
    DistributeableMeshObject,
    topoSetList
> TopoSetList;


/*---------------------------------------------------------------------------*\
                       Class topoSetList Declaration
\*---------------------------------------------------------------------------*/

class topoSetList
:
    public TopoSetList
{
public:
    enum SelectionType
    {
        ALL,
        INTERNAL,
        INTERFACE,
        BOUNDARY,
        INTERFACE_AND_BOUNDARY
    };
    static const NamedEnum<SelectionType, 5> selectionNames;

private:
// Private data

    //- Zones and sets are all stored as topo sets so they can be used
    //  interchangeably
    HashPtrTable<topoSet> cellTopoSets_;
    HashPtrTable<topoSet> faceTopoSets_;
    HashPtrTable<topoSet> pointTopoSets_;

    //- Hashed sets of toposet types
    wordHashSet cellSets_;
    wordHashSet cellZones_;

    wordHashSet faceSets_;
    wordHashSet faceZones_;

    wordHashSet pointSets_;
    wordHashSet pointZones_;


public:

    //- TypeName
    TypeName("topoSetList");

    //- Constructor
    topoSetList(const fvMesh& mesh);

    //- Destructor
    virtual ~topoSetList();

    // Public member functions

        //- Return cell sets
        const HashPtrTable<topoSet>& cellSets() const
        {
            return cellTopoSets_;
        }

        //- Return face sets
        const HashPtrTable<topoSet>& faceSets() const
        {
            return faceTopoSets_;
        }

        //- Return point sets
        const HashPtrTable<topoSet>& pointSets() const
        {
            return pointTopoSets_;
        }

        //- Update cell sets
        void updateCells
        (
            const dictionary& dict,
            const labelList& selectedCells,
            const bool isZone = false
        );

        //- Update face sets
        void updateFaces
        (
            const dictionary& dict,
            const labelList& selectedFaces,
            const boolList& flipMap,
            const bool isZone = false
        );

        //- Update points sets
        void updatePoints
        (
            const dictionary& dict,
            const labelList& selectedPoints,
            const bool isZone = false
        );

        //- Update sets
        void update
        (
            const dictionary& dict,
            const labelList& selectedCells,
            const labelList& selectedFaces,
            const boolList& flipMap,
            const labelList& selectedPoints
        );

        //- Modify the given sets
        void modifyTopoSet
        (
            HashPtrTable<topoSet>& topoSets,
            wordHashSet& zones,
            wordHashSet& sets,
            const topoSet& selected,
            const dictionary& dict
        );

        //- Extract interface cells
        labelList extractInterfaceCells(const labelList& cells) const;

        //- Extract interface cells
        autoPtr<topoSet> extractInterfaceCells(const topoSet& cells) const;

        //- Extract selected faces (return indices to keep)
        labelList extractSelectedFaces
        (
            const dictionary& dict,
            const labelList& faces,
            const bool defaultAll = false
        ) const;

        //- Extract selected faces after getting the kept indices
        //  This is important because the flipMap needs to be kept if a zone
        //  is used
        autoPtr<topoSet> extractSelectedFaces
        (
            const dictionary& dict,
            const topoSet& faces,
            const bool defaultAll = false
        ) const;

        //- Extract selected points
        labelList extractSelectedPoints
        (
            const dictionary& dict,
            const labelList& points,
            const bool defaultAll = false
        ) const;

        //- Extract selected points
        autoPtr<topoSet> extractSelectedPoints
        (
            const dictionary& dict,
            const topoSet& points,
            const bool defaultAll = false
        ) const;

        //- Clear all sets
        void clear();

        virtual bool movePoints()
        {
            return false;
        }

        virtual void reorderPatches
        (
            const labelUList& newToOld,
            const bool validBoundary
        );

        virtual void addPatch(const label patchi)
        {}

        //- Update any stored data for new labels. Not implemented.
        virtual void updateMesh(const mapPolyMesh& morphMap);

        //- Distribute sets
        virtual void distribute(const mapDistributePolyMesh& map);

        //- Send the cellZones to the mesh, by default zones are removed here
        void transferZones(const bool remove = true);

    // Write sets and zones
    bool writeSets() const;


    // Static functions

        //- Extract interface cells
        static labelList extractInterfaceCells
        (
            const polyMesh&,
            const labelList&
        );

        //- Extract interface cells
        static autoPtr<topoSet> extractInterfaceCells
        (
            const polyMesh& mesh,
            const topoSet& cells
        );

        //- Extract selected faces (return indices to keep)
        static labelList extractSelectedFaces
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const labelList& faces,
            const bool defaultAll = false
        );

        //- Extract selected faces after getting the kept indices
        //  This is important because the flipMap needs to be kept if a zone
        //  is used
        static autoPtr<topoSet> extractSelectedFaces
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const topoSet& faces,
            const bool defaultAll = false
        );

        //- Extract selected points
        static labelList extractSelectedPoints
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const labelList& points,
            const bool defaultAll = false
        );

        //- Extract selected points
        static autoPtr<topoSet> extractSelectedPoints
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const topoSet& points,
            const bool defaultAll = false
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

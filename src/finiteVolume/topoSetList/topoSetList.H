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
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class topoSetList;
typedef MeshObject
<
    fvMesh,
    UpdateableMeshObject,
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
            const labelList& selectedPoints
        );

        //- Modify the given sets
        void modifyTopoSet
        (
            const word& setType,
            const bool isZone,
            HashPtrTable<topoSet>& topoSets,
            wordHashSet& zones,
            wordHashSet& sets,
            const labelList& selected,
            const dictionary& dict
        );

        //- Extract interface cells
        labelList extractInterfaceCells(const labelList& cells) const;

        //- Extract selected faces
        labelList extractSelectedFaces
        (
            const dictionary& dict,
            const labelList& faces,
            const bool defaultAll = false
        ) const;

        //- Extract selected points
        labelList extractSelectedPoints
        (
            const dictionary& dict,
            const labelList& points,
            const bool defaultAll = false
        ) const;

        //- Clear all sets
        void clear();

        //- Update any stored data for new labels. Not implemented.
        virtual void updateMesh(const mapPolyMesh& morphMap);

        //- Move points
        virtual bool movePoints()
        {
            return false;
        }

        //- Send the cellZones to the mesh, by default zones are removed here
        void transferZones(const bool remove = true);

    // Write sets and zones
    bool writeSets() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::hexRefLocationMapper

Description
    Various for reading/decomposing/reconstructing/distributing refinement
    data.

SourceFiles
    hexRefLocationMapper.C

\*---------------------------------------------------------------------------*/

#ifndef hexRefLocationMapper_H
#define hexRefLocationMapper_H

#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "autoPtr.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class hexRefLocationMapper Declaration
\*---------------------------------------------------------------------------*/

class hexRefLocationMapper
{

private:

    // Private data

        //- Reference to the mesh
        const polyMesh& mesh_;

        //- Does this actually need to be used
        bool constructMap_;

        //- Saved data from the old mesh
        autoPtr<pointField> pointsOldPtr_;
        autoPtr<edgeList> edgesOldPtr_;
        autoPtr<faceList> facesOldPtr_;
        autoPtr<cellList> cellsOldPtr_;

        enum splitType
        {
            EDGE,
            FACE,
            CELL
        };

        DynamicList<label> edgeSplits_;
        DynamicList<label> newEdgeIndices_;

        DynamicList<label> faceSplits_;
        DynamicList<label> newFaceIndices_;

        DynamicList<label> cellSplits_;
        DynamicList<label> newCellIndices_;

        //- Interpolation weights
        mutable autoPtr<scalarListList> edgeWeightsPtr_;
        mutable autoPtr<scalarListList> faceWeightsPtr_;
        mutable autoPtr<scalarListList> cellWeightsPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        hexRefLocationMapper(const hexRefLocationMapper&);

        //- Disallow default bitwise assignment
        void operator=(const hexRefLocationMapper&);

        const scalarListList& edgeWeights() const;
        const scalarListList& faceWeights() const;
        const scalarListList& cellWeights() const;


public:

    //- Construct from mesh
    hexRefLocationMapper(const polyMesh& mesh);

    //- Destructor
    ~hexRefLocationMapper();


    // Member Functions

        //- Set if the map needs to be constructed
        void needMap(const bool need = true)
        {
            constructMap_ = need;
        }

        //- Add split edges
        void addSplitEdges(const labelList&, const labelList&);

        //- Add split faces
        void addSplitFaces(const labelList&, const labelList&);

        //- Add split cells
        void addSplitCells(const labelList&, const labelList&);

        //- Correct the given point field
        template<class Type>
        void interpolate(Field<Type>&) const;

        //- Interpolate to the given edge/face/cell centres
        void interpolateMidPoints(pointField&) const;

        //- Clear the mapping
        void clearOut();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "hexRefLocationMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

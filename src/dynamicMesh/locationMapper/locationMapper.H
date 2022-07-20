/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applies Technologies
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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
    Foam::locationMapper

Description
    Class for handeling the mapping of locations after refinement. This must
    be modified by the mesh refiner due to the fact that determining the group
    of averaged points that new face centre

SourceFiles
    locationMapper.C
    locationMapperTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef locationMapper_H
#define locationMapper_H

#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "autoPtr.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class locationMapper Declaration
\*---------------------------------------------------------------------------*/

class locationMapper
{

private:

    // Private data

        //- Reference to the mesh
        const polyMesh& mesh_;

        //- Does this actually need to be used
        mutable bool constructMap_;

        //- Saved data from the old mesh
        mutable autoPtr<pointField> pointsOldPtr_;
        mutable autoPtr<edgeList> edgesOldPtr_;
        mutable autoPtr<faceList> facesOldPtr_;
        mutable autoPtr<cellList> cellsOldPtr_;

        enum splitType
        {
            EDGE,
            FACE,
            CELL
        };

        mutable DynamicList<label> edgeSplits_;
        mutable DynamicList<label> newEdgeIndices_;

        mutable DynamicList<label> faceSplits_;
        mutable DynamicList<label> newFaceIndices_;

        mutable DynamicList<label> cellSplits_;
        mutable DynamicList<label> newCellIndices_;

        //- Interpolation weights
        mutable autoPtr<scalarListList> edgeWeightsPtr_;
        mutable autoPtr<scalarListList> faceWeightsPtr_;
        mutable autoPtr<scalarListList> cellWeightsPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        locationMapper(const locationMapper&);

        //- Disallow default bitwise assignment
        void operator=(const locationMapper&);

        const scalarListList& edgeWeights() const;
        const scalarListList& faceWeights() const;
        const scalarListList& cellWeights() const;


public:

    //- Construct from mesh
    locationMapper(const polyMesh& mesh);

    //- Destructor
    ~locationMapper();


    // Member Functions

        //- Set if the map needs to be constructed
        void needMap(const bool need = true) const
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
        void clearOut() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "locationMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::mappedPatchSelector

Description
    Class used to return mapping functions using either a mappedPolyPatch
    or a mappedMovingPolyPatch


Note

SourceFiles
    mappedPatchSelector.C

\*---------------------------------------------------------------------------*/

#ifndef mappedPatchSelector_H
#define mappedPatchSelector_H

#include "fvPatch.H"
#include "mappedPatchBase.H"
#include "mappedMovingPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyPatch;
class polyMesh;
class mapDistribute;

/*---------------------------------------------------------------------------*\
                       Class mappedPatchSelector Declaration
\*---------------------------------------------------------------------------*/

class mappedPatchSelector
{
protected:

    // Protected data

        //- Pointer to mappedPatchBase
        const mappedPatchBase* mappedPatchPtr_;

        //- Pointer to mappedMovingPatchBase
        const mappedMovingPatchBase* mappedMovingPatchPtr_;


public:

    //- Runtime type information
    TypeName("mappedPatchSelector");


    // Constructors

        //- Construct from patch
        mappedPatchSelector(const fvPatch&);


    //- Destructor
    virtual ~mappedPatchSelector();


    // Member Functions

        // Access

            //- Region to sample
            inline const word& sampleRegion() const;

            //- Patch (only if NEARESTPATCHFACE)
            inline const word& samplePatch() const;

            //- PatchGroup (only if NEARESTPATCHFACE)
            inline const word& coupleGroup() const;

            //- Return size of mappedMoving mesh/patch/boundary
            inline label sampleSize() const;

            //- Cached sampleRegion != mesh.name()
            inline bool sameRegion() const;

            //- Return reference to the parallel distribution map
            inline const mapDistribute& map() const;

            //- Return reference to the AMI interpolator
            inline const AMIInterpolation& AMI
            (
                const bool forceUpdate = false
            ) const;

            //- Return a pointer to the AMI projection surface
            const autoPtr<Foam::searchableSurface>& surfPtr() const;

            //- Get the region mesh
            const polyMesh& sampleMesh() const;

            //- Get the patch on the region
            const polyPatch& samplePolyPatch() const;


        // Helpers

            //- Get the sample points
            tmp<pointField> samplePoints() const;

            //- Get a point on the face given a face decomposition method:
            //  face-centre-tet : face centre. Returns index of face.
            //  face-planes     : face centre. Returns index of face.
            //  face-diagonal   : intersection of ray from cellcentre to
            //                    facecentre with any of the triangles.
            //                    Returns index (0..size-2) of triangle.
            pointIndexHit facePoint
            (
                const polyMesh&,
                const label facei,
                const polyMesh::cellDecomposition
            );


        // Distribute

            //- Wrapper around map/interpolate data distribution
            template<class Type>
            void distribute(List<Type>& lst) const;

            //- Wrapper around map/interpolate data distribution with operation
            template<class Type, class CombineOp>
            void distribute(List<Type>& lst, const CombineOp& cop) const;

            //- Wrapper around map/interpolate data distribution
            template<class Type>
            void reverseDistribute(List<Type>& lst) const;

            //- Wrapper around map/interpolate data distribution with operation
            template<class Type, class CombineOp>
            void reverseDistribute(List<Type>& lst, const CombineOp& cop) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mappedPatchSelectorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedPatchSelectorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

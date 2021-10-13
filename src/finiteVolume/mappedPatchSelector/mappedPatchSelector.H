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
    Foam::mappedPatchSelector

Description
    Class used to return mapping functions using either a
    mappedPolyPatch or a mappedMovingPolyPatch


Note

SourceFiles
    mappedPatchSelector.C

\*---------------------------------------------------------------------------*/

#ifndef mappedPatchSelector_H
#define mappedPatchSelector_H

#include "polyPatch.H"
#include "mappedPatchBase.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class mapDistribute;

/*---------------------------------------------------------------------------*\
                       Class mappedPatchSelector Declaration
\*---------------------------------------------------------------------------*/

class mappedPatchSelector
{
protected:

    // Protected data

        //- Reference to the poly patch
        const polyPatch& patch_;

        //- Pointer to mappedPatchBase
        const mappedPatchBase* mappedPatchPtr_;

        //- Interpolation from face to point values
        mutable autoPtr<primitivePatchInterpolation> pointInterpolatorPtr_;


public:

    //- Runtime type information
    TypeName("mappedPatchSelector");


    // Constructors

        //- Construct from patch
        mappedPatchSelector(const polyPatch&);


    //- Destructor
    virtual ~mappedPatchSelector();


    // Member Functions

        // Access

            //- Get this mesh
            const polyMesh& mesh() const;

            //- Get the this patch
            const polyPatch& patch() const;

            //- Get the region mesh
            const polyMesh& sampleMesh() const;

            //- Get the patch on the region
            const polyPatch& samplePolyPatch() const;

            //- Return the point interpolator
            const primitivePatchInterpolation& pointInterpolator() const;


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


        // Point interpolation

            //- Wrapper around map/interpolate data distribution
            template<class Type>
            void distributePoint(List<Type>& lst) const;

            //- Wrapper around map/interpolate data distribution with operation
            template<class Type, class CombineOp>
            void distributePoint
            (
                List<Type>& lst,
                const CombineOp& cop
            ) const;


        // Clear the weights
        void clearOut();
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

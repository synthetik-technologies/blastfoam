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
    Foam::wallCollisionData

Description

SourceFiles
    wallCollisionData.C

\*---------------------------------------------------------------------------*/

#ifndef immersedwallCollisionData_H
#define immersedwallCollisionData_H

#include "labelPair.H"
#include "immersedBoundaryObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                     Class wallCollisionData Declaration
\*---------------------------------------------------------------------------*/

class wallCollisionData
:
    public labelPairList
{
    // Private Data

        //- Wall patch
        const polyPatch& patch_;

        //- Patch bounding box (modified outside)
        boundBox patchBB_;

        //- Immersed boundary object
        const immersedBoundaryObject& object_;

        //- Weights
        scalarField weights_;

        //- Reduced weights
        scalarField reducedWeights_;

        //- Average hit point
        vector hitPoint_;

        //- Average normal from wall
        vector normal_;

        //- Average relative velocity
        vector v_;

        //- Saved velocity of the patch
        vectorField vW_;

        //- Face of the immersedBoundaryObject closest to the hit point
        //  (smallest distance from the centre of rotation)
        //  In parallel this is -1 if not on the proc
        label hitPointIndex_;

        //- Reduced list of label pair that has had positive relative
        //  velocities removed
        labelPairList reducedMap_;


public:

    // Constructors

        //- Construct from components
        wallCollisionData
        (
            const polyPatch& patch,
            const immersedBoundaryObject& object
        );

    //- Destructor
    virtual ~wallCollisionData();


    // Member Functions

         //- Calculate the weights based on overlap distance
        //  for object-object relations and update the map to only include
        // points with negative relative normal velocities
        void calcMapping();

        //- Return the weights
        const scalarField& weights() const
        {
            return weights_;
        }

        //- Return the reduced weights
        const scalarField& reducedWeights() const
        {
            return reducedWeights_;
        }

        //- Return the hit point
        const vector& hitPoint() const
        {
            return hitPoint_;
        }

        //- Return the surface normal
        const vector& normal() const
        {
            return normal_;
        }

        //- Return the relative velocity
        const vector& v() const
        {
            return v_;
        }

        //- Return the hit index
        label hitPointIndex() const
        {
            return hitPointIndex_;
        }

        //- Return the reduced map
        const labelPairList& reducedMap() const
        {
            return reducedMap_;
        }

        //- Return the patch
        const polyPatch& patch() const
        {
            return patch_;
        }

        //- Return the patch bounding box
        const boundBox& patchBB() const
        {
            return patchBB_;
        }

        //- Access the patch bounding box
        boundBox& patchBB()
        {
            return patchBB_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

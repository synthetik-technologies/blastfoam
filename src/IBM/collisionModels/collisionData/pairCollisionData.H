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
    Foam::pairCollisionData

Description
    Class for calculating and storing pair collision data.

SourceFiles
    pairCollisionData.C

\*---------------------------------------------------------------------------*/

#ifndef immersedpairCollisionData_H
#define immersedpairCollisionData_H

#include "labelPair.H"
#include "immersedBoundaryObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                     Class pairCollisionData Declaration
\*---------------------------------------------------------------------------*/

class pairCollisionData
:
    public labelPairList
{
    // Private Data

        //- First immersedBoundaryObject
        const immersedBoundaryObject& objectA_;

        //- Second immersedBoundaryObject
        const immersedBoundaryObject& objectB_;

        //- Weights
        scalarField weights_;

        //- Reduced weights, i.e. positive normal velocities have been removed
        scalarField reducedWeights_;

        //- Average hit point
        vector hitPoint_;

        //- Average normal from objectA
        vector normalA_;

        //- Average normal from objectB
        vector normalB_;

        //- Average relative velocity (vA - vB)
        vector vAB_;

        //- Face of the immersedBoundaryObject closest to the hit point
        //  (smallest distance from the centre of rotation)
        //  In parallel this is -1 if not on the proc
        label hitPointIndexA_;
        label hitPointIndexB_;

        //- Reduced list of label pair that has had positive relative
        //  velocities removed
        labelPairList reducedMap_;


public:

    // Constructors

        //- Construct from components
        pairCollisionData
        (
            const immersedBoundaryObject& OA,
            const immersedBoundaryObject& OB
        );

    //- Destructor
    virtual ~pairCollisionData();


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

        //- Return the surface normal of object A
        const vector& normalA() const
        {
            return normalA_;
        }

        //- Return the surface normal of object B
        const vector& normalB() const
        {
            return normalB_;
        }

        //- Return the relative velocity
        const vector& vAB() const
        {
            return vAB_;
        }

        //- Return the hit index
        label hitPointIndexA() const
        {
            return hitPointIndexA_;
        }

        //- Return the hit index
        label hitPointIndexB() const
        {
            return hitPointIndexB_;
        }

        //- Return the reduced map
        const labelPairList& reducedMap() const
        {
            return reducedMap_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

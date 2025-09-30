/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::localLeastSquaresVectors

Description
    local molecule least-squares gradient scheme vectors

SourceFiles
    localLeastSquaresVectors.C

\*---------------------------------------------------------------------------*/

#ifndef localLeastSquaresVectors_H
#define localLeastSquaresVectors_H

#include "MeshObjects.H"
#include "DemandDrivenMeshObject.H"
#include "fvMesh.H"
#include "surfaceFieldsFwd.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                     Class localLeastSquaresVectors Declaration
\*---------------------------------------------------------------------------*/

class localLeastSquaresVectors
:
    public DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        localLeastSquaresVectors
    >
{
    // Private data

        //- Use purely local gradient
        bool local_;

        //- Inverse of distance matrix
        mutable tensorField* AinvPtr_;

        //- Inverse of local distance matrix
        mutable tensorField* AinvLocalPtr_;


    // Private member functions

        //- Construct Least-squares gradient vectors
        void makeLeastSquaresVectors() const;

        //- Construct Least-squares gradient vectors
        void makeLocalLeastSquaresVectors() const;

        //- Stabilised  tensor inverse
        static tensor stabInv(const tensor&);


public:

    // Declare name of the class and its debug switch
    TypeName("localLeastSquaresVectors");


    // Constructors

        //- Construct given an fvMesh and the minimum determinant criterion
        localLeastSquaresVectors
        (
            const fvMesh&,
            const bool local
        );


    // Destructor

        virtual ~localLeastSquaresVectors();


    // Member functions

        //- Return reference to owner least square vectors
        const tensorField& Ainv() const;

        //- Return reference to neighbour least square vectors
        const tensorField& AinvLocal() const;


        //- Delete the least square vectors when the mesh moves
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

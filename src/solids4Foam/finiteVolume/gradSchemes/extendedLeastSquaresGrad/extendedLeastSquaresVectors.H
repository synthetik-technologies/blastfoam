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
    Foam::extendedLeastSquaresVectors

Description
    Extended molecule least-squares gradient scheme vectors

SourceFiles
    extendedLeastSquaresVectors.C

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAMESIORFOUNDATION

#ifndef extendedLeastSquaresVectors_H
#define extendedLeastSquaresVectors_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "surfaceFieldsFwd.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class extendedLeastSquaresVectors Declaration
\*---------------------------------------------------------------------------*/

class extendedLeastSquaresVectors
:
    public MeshObject<fvMesh, MoveableMeshObject, extendedLeastSquaresVectors>
{
    // Private data

        //- Minimum determinant criterion to choose extra cells
        scalar minDet_;

        //- Least-squares gradient vectors
        mutable surfaceVectorField* pVectorsPtr_;
        mutable surfaceVectorField* nVectorsPtr_;

        mutable List<labelPair>* additionalCellsPtr_;
        mutable vectorField* additionalVectorsPtr_;


    // Private member functions

        //- Construct Least-squares gradient vectors
        void makeLeastSquaresVectors() const;


public:

    // Declare name of the class and its debug switch
    TypeName("extendedLeastSquaresVectors");


    // Constructors

        //- Construct given an fvMesh and the minimum determinant criterion
        extendedLeastSquaresVectors
        (
            const fvMesh&,
            const scalar minDet
        );


    // Destructor

        virtual ~extendedLeastSquaresVectors();


    // Member functions

        //- Return reference to owner least square vectors
        const surfaceVectorField& pVectors() const;

        //- Return reference to neighbour least square vectors
        const surfaceVectorField& nVectors() const;


        //- Return reference to additional cells
        const List<labelPair>& additionalCells() const;

        //- Return reference to additional least square vectors
        const vectorField& additionalVectors() const;


        //- Delete the least square vectors when the mesh moves
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

#endif // end of #ifdef OPENFOAMESIORFOUNDATION

// ************************************************************************* //

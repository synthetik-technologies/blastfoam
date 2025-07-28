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
    Foam::fv::localLeastSquaresGrad

Description
    Second-order gradient scheme using least-squares.

SourceFiles
    localLeastSquaresGrad.C

\*---------------------------------------------------------------------------*/

#ifndef localLeastSquaresGrad_H
#define localLeastSquaresGrad_H

#include "gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class localLeastSquaresGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class localLeastSquaresGrad
:
    public fv::gradScheme<Type>
{
    // Private Data

        //- Minimum determinant criterion to choose extra cells
        bool local_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        localLeastSquaresGrad(const localLeastSquaresGrad&);

        //- Disallow default bitwise assignment
        void operator=(const localLeastSquaresGrad&);


public:

    //- Runtime type information
    TypeName("localLeastSquares");


    // Constructors

        //- Construct from Istream
        localLeastSquaresGrad(const fvMesh& mesh, Istream& schemeData)
        :
            gradScheme<Type>(mesh),
            local_(Switch(schemeData))
        {}


    // Member Functions

        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "localLeastSquaresGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

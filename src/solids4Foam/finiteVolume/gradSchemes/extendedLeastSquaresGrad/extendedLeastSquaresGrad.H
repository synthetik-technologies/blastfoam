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
    Foam::fv::extendedLeastSquaresGrad

Description
    Second-order gradient scheme using least-squares.

SourceFiles
    extendedLeastSquaresGrad.C

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAMESIORFOUNDATION

#ifndef extendedLeastSquaresGrad_H
#define extendedLeastSquaresGrad_H

#include "gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class extendedLeastSquaresGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class extendedLeastSquaresGrad
:
    public fv::gradScheme<Type>
{
    // Private Data

        //- Minimum determinant criterion to choose extra cells
        scalar minDet_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        extendedLeastSquaresGrad(const extendedLeastSquaresGrad&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const extendedLeastSquaresGrad&) = delete;


public:

    //- Runtime type information
    TypeName("extendedLeastSquares");


    // Constructors

        //- Construct from Istream
        extendedLeastSquaresGrad(const fvMesh& mesh, Istream& schemeData)
        :
            gradScheme<Type>(mesh),
            minDet_(readScalar(schemeData))
        {
            if (minDet_ < 0 || minDet_ > 8)
            {
                FatalIOErrorIn
                (
                    "extendedLeastSquaresGrad"
                    "(const fvMesh&, Istream& schemeData)",
                    schemeData
                )   << "Minimum determinant = " << minDet_
                    << " should be >= 0 and <= 8"
                    << exit(FatalIOError);
            }
        }


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
#   include "extendedLeastSquaresGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

#endif // end of #ifdef OPENFOAMESIORFOUNDATION

// ************************************************************************* //

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
    RBFInterpolation

Description
    Radial basis function interpolation class

Description
    Interpolation class which uses Radial Basis Functions to interpolate
    field from given data points to arbitrary set of points.

    The coefficient vectors, alpha and beta are determined by solving
    the system:

    | db | = | Mbb Pb | | alpha |
    | 0  |   | Pb  0  | |  beta |

    where db are the given field values at data carrier points.
    Mbb the boundary RBF correlation matrix (NbxNb), containing RBF evaluations
    at the boundary nodes, and Pb some linear polynomial matrix (Nbx4).

    In cases where far field data is not of interest, a cutoff function
    is used to eliminate unnecessary data points in the far field

Author
    Frank Bos, TU Delft.  All rights reserved.
    Dubravko Matijasevic, FSB Zagreb.
    Reorganisation by Hrvoje Jasak, Wikki Ltd.

    David Blom, TU Delft. All rights reserved.

SourceFiles
    RBFInterpolation.C
    RBFInterpolationTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef RBFInterpolation_H
#define RBFInterpolation_H

#include "patchToPatchMapping.H"
#include "dictionary.H"
#include "RBFFunction.H"
#include "point.H"


#ifdef USE_EIGEN

#ifndef FULLDEBUG
#define NDEBUG
#endif

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/QR>

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class RBFInterpolation Declaration
\*---------------------------------------------------------------------------*/

class RBFInterpolation
{
#ifdef USE_EIGEN
    typedef Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::FullPivLU<Matrix> LUMatrix;
    typedef Eigen::LLT<Matrix> LLTMatrix;
    typedef Eigen::ColPivHouseholderQR<Matrix> QRMatrix;

    template<class Type>
    using VectorMatrix =
        Eigen::Map
        <
            Eigen::Matrix
            <
                scalar,
                Eigen::Dynamic,
                Eigen::Dynamic,
                Eigen::RowMajor
            >
        >;

    template<class Type>
    using CVectorMatrix =
        Eigen::Map
        <
            const Eigen::Matrix
            <
                scalar,
                Eigen::Dynamic,
                Eigen::Dynamic,
                Eigen::RowMajor
            >
        >;

    template<class Type>
    using Block = Eigen::Block<Type>;
#endif

    // Private data

        patchToPatchMapping::Mapping mapping_;

        //- Reference to control points
        const vectorField& controlPoints_;

        //- Reference to all points
        const vectorField& dataPoints_;

        //- RBF function
        autoPtr<RBFFunction> RBF_;

#ifdef USE_EIGEN
        //- Positions using eigen classes
        CVectorMatrix<vector> positions_;
        CVectorMatrix<vector> positionsInterpolation_;

        //- Evaluation matrix (output x input)
        mutable Matrix* APtr_;

        //- Decomposition of the interpolation matrix
        mutable LLTMatrix* decompLLTPtr_;
        mutable QRMatrix* decompQRPtr_;

        //- Matrix used for conservative interpolation
        mutable Matrix* HhatPtr_;

        //- Polynomial matrix of the input mesh (for separate polynomial)
        mutable Matrix* QPtr_;

        //- Polynomial matrix of the output mesh (for separate polynomial)
        mutable Matrix* VPtr_;

        //- Decomposition of the polynomial (for separate polynomial)
        mutable QRMatrix* QRPtr_;
#endif

        bool separate_;

        //- Add polynomials to RBF matrix
        bool polynomials_;

        //- Dead directions
        Vector<label> validDirs_;
        label nDead_;

        //- Polynomial size
        label polySize_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const RBFInterpolation&);

#ifdef USE_EIGEN
        //- Fill polynomial entries in the given matrix
        void fillPolynomialEntries
        (
            Matrix& m,
            const label size,
            const UList<vector>& pts
        ) const;


        //- Calculate evaluation matrix
        void calcA() const;

        //- Return evaluation matrix
        const Matrix& A() const;

        //- Calculate interpolation matrices
        void calc() const;

        //- Return decomposition of the interpolation matrix
        const LLTMatrix& decompLLT() const;
        const QRMatrix& decompQR() const;

        //- Conservative interpolation matrix
        const Matrix& Hhat() const;

        //- Calculate conservative interpolation matrix
        void calcHhat() const;

        //- Polynomial matrix of the input mesh (for separate polynomial)
        const Matrix& Q() const;

        //- Polynomial matrix of the output mesh (for separate polynomial)
        const Matrix& V() const;

        //- Decomposition of the polynomial (for separate polynomial)
        const QRMatrix& QR() const;
#endif

        //- Return pointer to data
        template<class T>
        static const scalar* getData(const UList<T>& lst)
        {
            return lst.cdata()->v_;
        }

        //- Access pointer to data
        template<class T>
        static scalar* getData(UList<T>& lst)
        {
            return lst.data()->v_;
        }

        //- Clear out
        void clearOut();


public:

    ClassName("RBFInterpolation");

    // Constructors

        //- Construct from components
        RBFInterpolation
        (
            const dictionary& dict,
            const vectorField& controlPoints,
            const vectorField& dataPoints,
            const patchToPatchMapping::Mapping mapping
        );

        //- Construct from components
        RBFInterpolation
        (
            const word& type,
            const dictionary& dict,
            const vectorField& controlPoints,
            const vectorField& dataPoints,
            const patchToPatchMapping::Mapping mapping
        );

        //- Construct as copy
        RBFInterpolation(const RBFInterpolation&);


    // Destructor

        virtual ~RBFInterpolation();


    // Member Functions

        //- Return reference to control points
        const vectorField& controlPoints() const
        {
            return controlPoints_;
        }

        //- Reference to all points
        const vectorField& dataPoints() const
        {
            return dataPoints_;
        }

        //- Interpolate
        template<class Type>
        tmp<Field<Type>> interpolate(const Field<Type>& fromField) const;

        // //- Interpolate
        // template<class Type>
        // tmp<Field<Type>> interpolate2(const Field<Type>& fromField) const;

        //- Interpolate
        template<class Type>
        void interpolate
        (
            const Field<Type>& fromField,
            Field<Type>& toField
        ) const;

        // //- Interpolate
        // template<class Type>
        // void interpolate2
        // (
        //     const Field<Type>& fromField,
        //     Field<Type>& toField
        // ) const;

        //- Move points
        void movePoints();
};

template<>
const scalar* RBFInterpolation::getData(const UList<scalar>& lst);

template<>
scalar* RBFInterpolation::getData(UList<scalar>& lst);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "RBFInterpolationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

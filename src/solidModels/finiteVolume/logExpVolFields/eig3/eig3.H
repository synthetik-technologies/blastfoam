/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    eig3

Description
    Eigen-decomposition for symmetric 3x3 real matrices.
    Public domain, copied from the public domain Java library JAMA.

SourceFiles
    eig3.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef _eig_h
#define _eig_h

#include "tensor.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class eig3 Declaration
\*---------------------------------------------------------------------------*/

class eig3
{
public:

    // Static member functions

        //- Symmetric matrix A => eigenvectors in columns of V, corresponding
        //  eigenvalues in d.

        static void tql2(tensor& V, vector& d, vector& e);

        static double hypot2(double x, double y);

        static void tred2(tensor& V, vector& d, vector& e);

        static void eigen_decomposition(const tensor& A, tensor& V, vector& d);

        static void eigen_decomposition
        (
            const symmTensor& sA, tensor& V, vector& d
        );

        static void eigen_decomposition
        (
            const symmTensor& sA, tensor& V, diagTensor& d
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2019 Alberto Passalacqua
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    Foam::eigenSolver

Description
    Eigenvalues and eigenvectors of a square real matrix

Acknowledgments
    This implementation derives almost completely from TNT/JAMA, a public-domain
    library developed at NIST, available at http://math.nist.gov/tnt/index.html
    Their implementation was based upon EISPACK.

    The tridiagonaliseSymmMatrix, symmTridiagQL, Hessemberg and realSchur
    methods are based on the Algol procedures tred2 by Bowdler, Martin, Reinsch,
    and Wilkinson, Handbook for Auto. Comp., Vol. II-Linear Algebra, and the
    corresponding FORTRAN subroutine in EISPACK.

SourceFiles
    eigenSolver.C

\*---------------------------------------------------------------------------*/

#ifndef eigenSolver_H
#define eigenSolver_H

#include "scalar.H"
#include "scalarMatrices.H"
#include "complex.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class eigenSolver Declaration
\*---------------------------------------------------------------------------*/

class eigenSolver
{
private:

    // Private data

        //- Real part of eigenvalues
        scalarDiagonalMatrix eigenvaluesRe_;

        //- Imaginary part of eigenvalues
        scalarDiagonalMatrix eigenvaluesIm_;

        //- Eigenvectors matrix
        scalarSquareMatrix eigenvectors_;

        //- Matrix with the same size of original matrix
        scalarSquareMatrix H_;

        //- Number of rows and columns
        label n_;

    // Private member functions

        //- Checks matrix for symmetry
        bool isSymmetric(const scalarSquareMatrix& A);

        //- Householder transform of a symmetric matrix to tri-diagonal form
        void tridiagonaliseSymmMatrix();

        //- Symmetric tri-diagonal QL algorithm
        void symmTridiagQL();

        //- Reduction of non-symmetric matrix to Hessemberg form
        void Hessenberg();

        //- Reduction from Hessenberg to real Schur form
        void realSchur();


public:

    // Constructors

        //- Construct from a scalarSquareMatrix
        eigenSolver(const scalarSquareMatrix& A);

        //- Construct from a scalarSquareMatrix and symmetry information
        //  Does not perform symmetric check
        eigenSolver(const scalarSquareMatrix& A, bool symmetric);

        //- Disallow default bitwise copy construct
        eigenSolver(const eigenSolver&) = delete;


    // Access Functions

    //- Return real part of the eigenvalues
        const scalarDiagonalMatrix& eigenvaluesRe() const
        {
            return eigenvaluesRe_;
        }

        const scalarDiagonalMatrix& eigenvaluesIm() const
        {
            return eigenvaluesIm_;
        }

        //- Return eigenvectors
        const scalarSquareMatrix& eigenvectors() const
        {
            return eigenvectors_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const eigenSolver&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

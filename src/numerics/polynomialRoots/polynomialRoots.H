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
    Solve the roots of a polynomial

Acknowledgments
    This implementation derives almost completely from TNT/JAMA, a public-domain
    library developed at NIST, available at http://math.nist.gov/tnt/index.html
    Their implementation was based upon EISPACK.

    The tridiagonaliseSymmMatrix, symmTridiagQL, Hessemberg and realSchur
    methods are based on the Algol procedures tred2 by Bowdler, Martin, Reinsch,
    and Wilkinson, Handbook for Auto. Comp., Vol. II-Linear Algebra, and the
    corresponding FORTRAN subroutine in EISPACK.

SourceFiles
    polynomialRoots.C

\*---------------------------------------------------------------------------*/

#ifndef polynomialRoots_H
#define polynomialRoots_H

#include "scalar.H"
#include "List.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class eigenSolver Declaration
\*---------------------------------------------------------------------------*/

class polynomialRoots
{
    // Private data

        // Polynomial coefficients
        List<scalar> coeffs_;

        // Saved real roots
        List<scalar> rootsRe_;

        // Saved imaginary roots
        List<scalar> rootsIm_;


public:

    // Constructors

        //- Construct from a scalarSquareMatrix
        polynomialRoots(const List<scalar>& P);

        //- Disallow default bitwise copy construct
        polynomialRoots(const polynomialRoots&) = delete;


    // Access Functions

        //- Return the polynomial coeffients
        const List<scalar>& coeffs() const
        {
            return coeffs_;
        }

        //- Return the function as a word
        word polyName() const;

        //- Evaluate the polynomial
        static scalar eval(const List<scalar>& coeffs, const scalar x);

        //- Evaluate the polynomial
        scalar eval(const scalar x) const;

        //- Return real parts of the roots
        const List<scalar>& rootsRe() const
        {
            return rootsRe_;
        }

        //- Return imaginary parts of the roots
        const List<scalar>& rootsIm() const
        {
            return rootsIm_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const polynomialRoots&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::polynomialRoots

Description
    Solve the roots of a polynomial

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
                    Class polynomialRoots Declaration
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
        static word polyName(const List<scalar>& coeffs);

        //- Return the function as a word
        word polyName() const;

        //- Evaluate the polynomial
        static scalar eval(const List<scalar>& coeffs, const scalar x);

        //- Evaluate the polynomial
        scalar eval(const scalar x) const;

        //- Evaluate the polynomial
        static List<scalar> derivative(const List<scalar>& coeffs);

        //- Evaluate the polynomial
        List<scalar> derivative() const;

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

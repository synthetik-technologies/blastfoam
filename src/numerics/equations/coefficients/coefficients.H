/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::coefficients

Description
    Base class for holding modifiable coefficents

\*---------------------------------------------------------------------------*/

#ifndef coefficients_H
#define coefficients_H

#include "List.H"
#include "scalar.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class coefficients Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class CoeffType>
class coefficients
{
protected:
    // Protected member data

        //- List of coefficients
        mutable Field<CoeffType> coeffs_;


public:

    // Constructors

        //- Construct given number of coeffients
        coefficients(const label nCoeffs);

        //- Construct given coeffients
        coefficients(const List<CoeffType>& coeffs);


    //- Destructor
    virtual ~coefficients();

    // Member functions

        //- Return the number of coefficients
        inline label n() const 
        {
            return coeffs_.size();
        }

        //- Return the coefficients
        const Field<CoeffType>& coeffs() const
        {
            return coeffs_;
        }

        //- Access the coefficients
        Field<CoeffType>& coeffsRef() const
        {
            return coeffs_;
        }

        virtual void coeffJ
        (
            const List<scalarList>& x,
            const label li,
            RectangularMatrix<Type>& J
        ) const
        {
            NotImplemented;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "coefficients.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

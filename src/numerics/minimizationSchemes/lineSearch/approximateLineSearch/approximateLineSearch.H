/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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
    Foam::approximateLineSearch

Description
    Approximate line search

SourceFiles
    approximateLineSearch.C

\*---------------------------------------------------------------------------*/

#ifndef approximateLineSearch_H
#define approximateLineSearch_H

#include "lineSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class approximateLineSearch Declaration
\*---------------------------------------------------------------------------*/

class approximateLineSearch
:
    public lineSearch
{
protected:
    // Private member data

        //- Initial step size
        scalar alpha0_;
        mutable scalar alpha_;

        //- Gradient control parameter (0-1)
        //  0: any change is acceptable
        //  1: Must decrease by more than first order approximation
        //  1e-4: default
        mutable scalar beta_;

        //- Step size increase factor
        scalar p_;


public:

    //- Runtime type information
    TypeName("approximate");


    // Constructors

        //- Construct for given equation
        approximateLineSearch
        (
            const scalarUnivariateEquation& eqns,
            const dictionary& dict
        );


    //- Destructor
    virtual ~approximateLineSearch();


    // Member Functions

        //- Line search using by limiting change in function value
        virtual void search
        (
            const scalarList& x0,
            const scalarList& grad,
            const label li,
            scalarList& xNew
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

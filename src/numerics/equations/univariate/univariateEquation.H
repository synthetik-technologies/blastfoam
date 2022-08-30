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
    Foam::univariateEquation

Description
    Abstract class for an equation that takes in a scalarList and outputs
    a single scalar value.

\*---------------------------------------------------------------------------*/

#ifndef univariateEquation_H
#define univariateEquation_H

#include "multivariateEquation.H"
#include "scalarField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class univariateEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class univariateEquation
:
    public multivariateEquation<Type>
{
public:

    typedef Type EqnType;

    // Constructors

        //- Null constructor
        univariateEquation()
        :
            multivariateEquation<Type>()
        {}

        //- Construct given name
        univariateEquation(const string& eqnString)
        :
            multivariateEquation<Type>(eqnString)
        {}

        //- Construct from dictionary
        univariateEquation(const dictionary& dict)
        :
            multivariateEquation<Type>(dict)
        {}


    //- Destructor
    virtual ~univariateEquation()
    {}


    //- Public member functions

        //- Return the function value
        virtual Type fX
        (
            const UList<scalar>& x,
            const label li
        ) const = 0;

        //- Return the function value
        virtual void dfdX
        (
            const UList<scalar>& x,
            const label li,
            List<Type>& dfdx
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

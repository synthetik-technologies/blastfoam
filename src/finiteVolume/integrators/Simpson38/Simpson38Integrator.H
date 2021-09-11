/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::Simpson38Integrator

Description
    Find root using the Simpson38 method

SourceFiles
    Simpson38Integrator.C

\*---------------------------------------------------------------------------*/

#ifndef Simpson38Integrator_H
#define Simpson38Integrator_H

#include "Integrator.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class Simpson38Integrator Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Simpson38Integrator
:
    public Integrator<Type>
{
// Private data

    const scalar oneThird = 1.0/3.0;
    const scalar twoThirds = 2.0/3.0;

public:

    //- Runtime type information
    TypeName("Simpson38");

    // Constructors

        //- Construct for given Equation
        Simpson38Integrator
        (
            const Equation<Type>& eqn,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        Simpson38Integrator(const Simpson38Integrator<Type>&) = delete;


    //- Destructor
    virtual ~Simpson38Integrator()
    {}


    // Member Functions

        //- Solve the root
        virtual Type integrate
        (
            const scalar x0,
            const scalar x1,
            const label li
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Simpson38Integrator<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "Simpson38Integrator.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

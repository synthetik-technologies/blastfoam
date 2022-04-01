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
    Foam::Simpson38Integrator

Description
    Simpson 3/8 integration method

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

    typedef univariateEquation<Type> equationType;

    const scalar oneThird = 1.0/3.0;
    const scalar twoThirds = 2.0/3.0;

        //- Solve the root
        virtual Type integrate_
        (
            const scalar dx,
            const Type& f0,
            const Type& f13,
            const Type& f23,
            const Type& f1
        ) const;

        //- Solve the root
        virtual Type integrate_
        (
            const Type& Q,
            const scalar x0,
            const scalar x1,
            const Type& f0,
            const Type& f13,
            const Type& f23,
            const Type& f1,
            const scalar tol,
            const label li
        ) const;

public:

    //- Runtime type information
    TypeName("Simpson38");

    // Constructors

        //- Construct for given Equation
        Simpson38Integrator
        (
            const equationType& eqn,
            const dictionary& dict
        );

        //- Construct for given Equation
        Simpson38Integrator
        (
            const equationType& eqn,
            const integrator& inter
        );

        //- Disallow default bitwise copy construction
        Simpson38Integrator(const Simpson38Integrator<Type>&) = delete;


    //- Destructor
    virtual ~Simpson38Integrator()
    {}

    // Member Functions

        //- Integrate the equation
        virtual Type integrate
        (
            const scalar x0,
            const scalar x1,
            const label li
        ) const;

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

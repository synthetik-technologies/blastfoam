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
    Foam::TrapezoidalIntegrator

Description
    Trapezoidal integration method

SourceFiles
    TrapezoidalIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef TrapezoidalIntegrator_H
#define TrapezoidalIntegrator_H

#include "Integrator.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class TrapezoidalIntegrator Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TrapezoidalIntegrator
:
    public Integrator<Type>
{
    typedef Equation<Type> equationType;

        //- Solve the root
        virtual Type integrate_
        (
            const scalar dx,
            const Type& f0,
            const Type& f1
        ) const;

        //- Solve the root
        virtual Type integrate_
        (
            const Type& Q,
            const scalar x0,
            const scalar x1,
            const Type& f0,
            const Type& f1,
            const scalar tol,
            const label li
        ) const;

public:

    //- Runtime type information
    TypeName("trapezoidal");

    // Constructors

        //- Construct for given Equation
        TrapezoidalIntegrator
        (
            const equationType& eqn,
            const dictionary& dict
        );

        //- Construct for given Equation
        TrapezoidalIntegrator
        (
            const equationType& eqn,
            const integrator& inter
        );

        //- Disallow default bitwise copy construction
        TrapezoidalIntegrator(const TrapezoidalIntegrator<Type>&) = delete;


    //- Destructor
    virtual ~TrapezoidalIntegrator()
    {}


    // Member Functions

        //- Integrate the equation
        virtual Type integrate
        (
            const scalar x0,
            const scalar x1,
            const label li
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const TrapezoidalIntegrator<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "TrapezoidalIntegrator.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

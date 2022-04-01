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
    Foam::BooleIntegrator

Description
    Simpson 3/8 integration method

SourceFiles
    BooleIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef BooleIntegrator_H
#define BooleIntegrator_H

#include "Integrator.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class BooleIntegrator Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BooleIntegrator
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
            const Type& f14,
            const Type& f12,
            const Type& f34,
            const Type& f1
        ) const;

        //- Solve the root
        virtual Type integrate_
        (
            const Type& Q,
            const scalar x0,
            const scalar x1,
            const Type& f0,
            const Type& f14,
            const Type& f12,
            const Type& f34,
            const Type& f1,
            const scalar tol,
            const label li
        ) const;

public:

    //- Runtime type information
    TypeName("Boole");

    // Constructors

        //- Construct for given Equation
        BooleIntegrator
        (
            const equationType& eqn,
            const dictionary& dict
        );

        //- Construct for given Equation
        BooleIntegrator
        (
            const equationType& eqn,
            const integrator& inter
        );

        //- Disallow default bitwise copy construction
        BooleIntegrator(const BooleIntegrator<Type>&) = delete;


    //- Destructor
    virtual ~BooleIntegrator()
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
        void operator=(const BooleIntegrator<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "BooleIntegrator.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

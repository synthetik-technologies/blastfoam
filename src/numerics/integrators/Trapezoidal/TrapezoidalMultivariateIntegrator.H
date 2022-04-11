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
    Foam::TrapezoidalMultivariateIntegrator

Description
    Mid-point integration

SourceFiles
    TrapezoidalMultivariateIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef TrapezoidalMultivariateIntegrator_H
#define TrapezoidalMultivariateIntegrator_H

#include "MultivariateIntegrator.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class TrapezoidalMultivariateIntegrator Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TrapezoidalMultivariateIntegrator
:
    public MultivariateIntegrator<Type>
{
    typedef univariateEquation<Type> equationType;

    virtual Type integrateFunc
    (
        const scalarList& x0,
        const scalarList& x1,
        const label li
    ) const;

    //- Recursive add to the sum
    virtual void addCorners
    (
        const label diri,
        const scalarList& dx,
        const label li,
        scalarList& x,
        Type& fx
    ) const;

public:

    //- Runtime type information
    TypeName("trapezoidal");

    // Constructors

        //- Construct for given Equation
        TrapezoidalMultivariateIntegrator
        (
            const equationType& eqn,
            const dictionary& dict
        );

        //- Construct for given Equation
        TrapezoidalMultivariateIntegrator
        (
            const equationType& eqn,
            const multivariateIntegrator& inter
        );

        //- Disallow default bitwise copy construction
        TrapezoidalMultivariateIntegrator
        (
            const TrapezoidalMultivariateIntegrator<Type>&
        ) = delete;


    //- Destructor
    virtual ~TrapezoidalMultivariateIntegrator()
    {}


    // Member Functions

        //- Integrate the equation
        virtual Type integrate
        (
            const scalarList& x0,
            const scalarList& x1,
            const label li
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const TrapezoidalMultivariateIntegrator<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "TrapezoidalMultivariateIntegrator.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

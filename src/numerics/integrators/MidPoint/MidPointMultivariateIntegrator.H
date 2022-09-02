/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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
    Foam::MidPointMultivariateIntegrator

Description
    Mid-point integration

SourceFiles
    MidPointMultivariateIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef MidPointMultivariateIntegrator_H
#define MidPointMultivariateIntegrator_H

#include "MultivariateIntegrator.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class MidPointMultivariateIntegrator Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MidPointMultivariateIntegrator
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


public:

    //- Runtime type information
    TypeName("midPoint");

    // Constructors

        //- Construct for given Equation
        MidPointMultivariateIntegrator
        (
            const equationType& eqn,
            const dictionary& dict
        );

        //- Construct for given Equation
        MidPointMultivariateIntegrator
        (
            const equationType& eqn,
            const multivariateIntegrator& inter
        );

        //- Disallow default bitwise copy construction
        MidPointMultivariateIntegrator
        (
            const MidPointMultivariateIntegrator<Type>&
        ) = delete;


    //- Destructor
    virtual ~MidPointMultivariateIntegrator()
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
        void operator=(const MidPointMultivariateIntegrator<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "MidPointMultivariateIntegrator.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

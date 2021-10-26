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
    Foam::MidPointIntegrator

Description
    Mid-point integration

SourceFiles
    MidPointIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef MidPointIntegrator_H
#define MidPointIntegrator_H

#include "Integrator.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class MidPointIntegrator Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MidPointIntegrator
:
    public Integrator<Type>
{
        typedef univariateEquation<Type> equationType;

public:

    //- Runtime type information
    TypeName("midPoint");

    // Constructors

        //- Construct for given Equation
        MidPointIntegrator
        (
            const equationType& eqn,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        MidPointIntegrator(const MidPointIntegrator<Type>&) = delete;


    //- Destructor
    virtual ~MidPointIntegrator()
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
        void operator=(const MidPointIntegrator<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "MidPointIntegrator.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

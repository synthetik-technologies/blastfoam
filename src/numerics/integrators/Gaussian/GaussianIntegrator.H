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
    Foam::GaussianIntegrator

Description
    Variable node Gaussian quadrature

SourceFiles
    GaussianIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef GaussianIntegrator_H
#define GaussianIntegrator_H

#include "Integrator.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class GaussianIntegrator Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class GaussianIntegrator
:
    public Integrator<Type>
{
// Private data

    typedef univariateEquation<Type> equationType;

    //- Gauss weights
    scalarList ws_;

    //- Normalized positions
    scalarList xs_;

    // Private member functions

        //- Set the weights and positions
        void setQuadrature(const label N);

        //- Solve the root
        virtual Type integrate_(const scalar dx, const PtrList<Type>& fs) const;

        //- Solve the root
        virtual Type integrate_
        (
            const Type& Q,
            const scalar x0,
            const scalar x1,
            const scalarList& xs,
            const PtrList<Type>& fs,
            const scalar tol,
            const label li
        ) const;

public:

    //- Runtime type information
    TypeName("Gaussian");

    // Constructors

        //- Construct for given Equation
        GaussianIntegrator
        (
            const equationType& eqn,
            const dictionary& dict
        );

        //- Construct for given Equation
        GaussianIntegrator
        (
            const equationType& eqn,
            const integrator& inter
        );

        //- Disallow default bitwise copy construction
        GaussianIntegrator(const GaussianIntegrator<Type>&) = delete;


    //- Destructor
    virtual ~GaussianIntegrator()
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
        void operator=(const GaussianIntegrator<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "GaussianIntegrator.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

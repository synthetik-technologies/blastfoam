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
    Foam::GaussianMultivariateIntegrator

Description
    Variable node Gaussian quadrature

SourceFiles
    GaussianMultivariateIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef GaussianMultivariateIntegrator_H
#define GaussianMultivariateIntegrator_H

#include "MultivariateIntegrator.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class GaussianMultivariateIntegrator Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class GaussianMultivariateIntegrator
:
    public MultivariateIntegrator<Type>
{
// Private data

    typedef univariateEquation<Type> equationType;

    //- Gauss weights
    List<scalarList> ws_;

    //- Normalized positions
    List<scalarList> xs_;


    // Private member functions

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
            const scalarList& x0,
            const scalarList& x1,
            const label li,
            labelList& gi,
            Type& fx
        ) const;

public:

    //- Runtime type information
    TypeName("Gaussian");

    // Constructors

        //- Construct for given Equation
        GaussianMultivariateIntegrator
        (
            const equationType& eqn,
            const dictionary& dict
        );

        //- Construct for given Equation
        GaussianMultivariateIntegrator
        (
            const equationType& eqn,
            const multivariateIntegrator& inter
        );

        //- Disallow default bitwise copy construction
        GaussianMultivariateIntegrator
        (
            const GaussianMultivariateIntegrator<Type>&
        ) = delete;


    //- Destructor
    virtual ~GaussianMultivariateIntegrator()
    {}

    // Member Functions

        //- Integrate the equation
        virtual Type integrate
        (
            const scalarList& x0,
            const scalarList& x1,
            const label li
        ) const;

        //- Disallow default bitwise assignment
        void operator=(const GaussianMultivariateIntegrator<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "GaussianMultivariateIntegrator.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

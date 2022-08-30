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
    Foam::Integrator

Description
    Generic class for integrating a given univariate equation.

SourceFiles
    multivariateIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateIntegrator_H
#define multivariateIntegrator_H

#include "integratorBase.H"
#include "Field.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class multivariateIntegrator Declaration
\*---------------------------------------------------------------------------*/

class multivariateIntegrator
:
    public integratorBase
{
protected:
    // Protected data

        //- Tolerance
        Field<scalar> tolerance_;

        //- Absolute tolerances
        Field<scalar> absTolerance_;

        //- Maximum number of domain splits
        List<label> maxSplits_;

        //- Number of intervals
        List<label> nIntervals_;

        mutable List<label> intervals_;
        mutable List<scalar> minDx_;

        //- Reset counters and minimum spacing
        void reset(const List<scalar>& dx) const;


public:

    //- Runtime type information
    TypeName("multivariateIntegrator");

    // Constructors

        //- Construct for from dictinoary
        multivariateIntegrator(const label n, const dictionary& dict);


    //- Destructor
    virtual ~multivariateIntegrator()
    {}


    // Member Functions

        //- Return if the integration is adaptive
        virtual bool adaptive() const = 0;

        //- Return the number of intervals
        List<label> nIntervals() const
        {
            return adaptive() ? intervals_ : nIntervals_;
        }

        //- Access the number of intervals to use
        virtual void setNIntervals(const List<label> nIntervals)
        {
            adaptive_ = false;
            nIntervals_ = nIntervals;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

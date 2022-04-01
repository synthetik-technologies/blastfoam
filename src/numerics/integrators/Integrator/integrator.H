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
    Foam::Integrator

Description
    Generic class for integrating a given univariate equation.

SourceFiles
    integrator.C

\*---------------------------------------------------------------------------*/

#ifndef integrator_H
#define integrator_H

#include "UnivariateEquation.H"
#include "typeInfo.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class integrator Declaration
\*---------------------------------------------------------------------------*/

class integrator
{
protected:
    // Protected data

        //- Is the integration adaptive
        bool adaptive_;

        //- Tolerance
        scalar tolerance_;

        //- Maximum number of domain splits
        label maxSplits_;

        //- Number of intervals
        label nIntervals_;

        mutable label evals_;
        mutable label intervals_;
        mutable scalar minDx_;

        //- Reset counters and minimum spacing
        void reset(const scalar dx) const;


public:

    //- Runtime type information
    TypeName("integrator");

    // Constructors

        //- Construct for from dictinoary
        integrator(const dictionary& dict);


    //- Destructor
    virtual ~integrator()
    {}


    // Member Functions

        //- Return if the integration is adaptive
        virtual bool adaptive() const = 0;

        //- Return the number of steps
        label nSteps() const
        {
            return evals_;
        }

        //- Return the number of intervals
        label nIntervals() const
        {
            return adaptive() ? intervals_ : nIntervals_;
        }

        //- Access the number of steps
        virtual void setAdaptive
        (
            const label maxSplits = 10,
            const scalar tol = 1e-6
        )
        {
            adaptive_ = true;
            maxSplits_ = maxSplits;
            tolerance_ = tol;
        }

        //- Access the number of intervals to use
        virtual void setNIntervals(const label nIntervals)
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

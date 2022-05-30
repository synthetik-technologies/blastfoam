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
    Foam::nonLinearLeastSquares

Description
    General class for finding equation coefficients of an equations

SourceFiles
    leastSquares.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinearLeastSquares_H
#define nonLinearLeastSquares_H

#include "coefficients.H"
#include "UnivariateEquationsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonLinearLeastSquares Declaration
\*---------------------------------------------------------------------------*/

class nonLinearLeastSquares
{
protected:
// Protected data

    //- Convergence tolerances
    scalar tolerance_;

    //- The maximum number of sub-steps allowed
    label maxSteps_;


public:

    //- Runtime type information
    TypeName("nonLinearLeastSquares");

    //- Constructor
    nonLinearLeastSquares
    (
        const scalar tol = 1e-6,
        const label maxSteps = 10
    );


    //- Destructor
    virtual ~nonLinearLeastSquares() 
    {}


    // Member Functions

        //- Find the root of the equation between xHigh and xLow
        virtual void findCoeffs
        (
            const scalarUnivariateEquation& eqns,
            const List<scalarList>& x,
            const scalarList& y,
            const label li = 0
        ) const;

        //- Find the root of the equation between xHigh and xLow
        void findCoeffs
        (
            const scalarUnivariateEquation& eqn,
            const List<scalarList>& x,
            const scalarList& y,
            const scalarList& w,
            const label li = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

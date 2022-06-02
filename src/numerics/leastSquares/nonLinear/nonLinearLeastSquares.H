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

#include "UnivariateEquationsFwd.H"
#include "coefficientEquationsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonLinearLeastSquares Declaration
\*---------------------------------------------------------------------------*/

class nonLinearLeastSquares
{
// Private data

    //- Convergence tolerances
    scalar tolerance_;

    //- The maximum number of sub-steps allowed
    label maxSteps_;

    //- Return a UList to the necessary data
    template<class VarType>
    static void setData(const VarType& x, UList<scalar>& xi);


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

        //- Update the coefficients of the given equations given
        //  a list of x and y values
        template<class VarType>
        void findCoeffs
        (
            scalarUnivariateEquation& eqn,
            const UList<VarType>& x,
            const scalarField& y,
            const label li = 0
        ) const;

        //- Update the coefficients of the given equations given
        //  a list of x and y values and associated weights
        template<class VarType>
        void findCoeffs
        (
            scalarUnivariateEquation& eqn,
            const UList<VarType>& x,
            const scalarField& y,
            const scalarField& w,
            const label li = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonLinearLeastSquaresTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

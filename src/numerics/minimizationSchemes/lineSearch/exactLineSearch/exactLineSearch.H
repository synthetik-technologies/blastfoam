/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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
    Foam::exactLineSearch

Description
    General class for minimizing multivariate, scalar equations

SourceFiles
    exactLineSearch.C

\*---------------------------------------------------------------------------*/

#ifndef exactLineSearch_H
#define exactLineSearch_H

#include "lineSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class univariateMinimizationScheme;

/*---------------------------------------------------------------------------*\
                    Class exactLineSearch Declaration
\*---------------------------------------------------------------------------*/

class exactLineSearch
:
    public lineSearch
{
protected:
    // Protected member data

        //- Line search minimization scheme
        mutable autoPtr<univariateMinimizationScheme> lineSearcher_;

        // Bracketing controls

            //- Maximum number of bracketing iterations
            label maxBracketIter_;

            //- Initial distance
            scalar dx0_;

            //- Step size increase factor
            scalar k_;


    // Protected member functions

            //- Bracket the line search
            void bracket
            (
                const label li,
                const scalar x0 = 0.0
            ) const;


public:

    //- Runtime type information
    TypeName("exact");


    // Constructors

        //- Construct for given exactLineSearch
        exactLineSearch
        (
            const scalarUnivariateEquation& eqns,
            const dictionary& dict
        );

        //- Copy constructor
        exactLineSearch
        (
            const scalarUnivariateEquation& eqns,
            const exactLineSearch& ls
        );

        //- Clone
        virtual autoPtr<lineSearch> clone
        (
            const scalarUnivariateEquation& eqns
        ) const
        {
            return autoPtr<lineSearch>(new exactLineSearch(eqns, *this));
        }


    //- Destructor
    virtual ~exactLineSearch();


    // Member Functions

        //- Line search using a univariate minimization scheme to
        //  find the minimum given the gradient
        virtual void search
        (
            const scalarList& x0,
            const scalarList& grad,
            const label li,
            scalarList& xNew
        ) const;

        //- Search the line for the minimum
        virtual void searchDir
        (
            const scalarList& x0,
            const scalarList& dir,
            const label li,
            scalarList& xNew
        ) const;


    // Static functions

        //- Bracket solution
        static void bracket
        (
            lineSearchEquation& eqn,
            const scalar x0,
            const label li,
            const scalar dx0,
            const scalar k,
            const label maxIter
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

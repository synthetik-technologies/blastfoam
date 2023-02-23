/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023
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
    Foam::lineSearch

Description
    Base class for line search methods

SourceFiles
    lineSearch.C

\*---------------------------------------------------------------------------*/

#ifndef lineSearch_H
#define lineSearch_H

#include "UnivariateEquationsFwd.H"
#include "lineSearchEquation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class lineSearch Declaration
\*---------------------------------------------------------------------------*/

class lineSearch
{
protected:

    // Protected data

        //- Reference to scalarEquation
        const scalarUnivariateEquation& eqns_;

        //- Line search equation
        mutable lineSearchEquation lsEqn_;


public:

    //- Runtime type information
    TypeName("lineSearch");

    // Declare run-time constructor selection tables
        declareRunTimeSelectionTable
        (
            autoPtr,
            lineSearch,
            dictionary,
            (const scalarUnivariateEquation& eqns, const dictionary& dict),
            (eqns, dict)
        );


    // Constructors

        //- Construct for given equation
        lineSearch
        (
            const scalarUnivariateEquation& eqns,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        lineSearch(const lineSearch&) = delete;


    // Selectors

        //- Select using a dictionary
        static autoPtr<lineSearch> New
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );


    //- Destructor
    virtual ~lineSearch();


    // Member Functions

        //- Search the line for the minimum
        virtual void search
        (
            const scalarList& x0,
            const scalarList& grad,
            const label li,
            scalarList& xNew
        ) const = 0;

        //- Return if the solution converged
        // virtual bool converged() const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const lineSearch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

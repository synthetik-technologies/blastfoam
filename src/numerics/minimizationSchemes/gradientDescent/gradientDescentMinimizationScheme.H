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
    Foam::gradientDescentMinimizationScheme

Description
    Minimization scheme using either gradient computed using a finite
    gradient of an implemented, analytical version

References:
    \verbatim
        "Algorithms for Optimization",
        Kochenderfer, M. J., Wheeler, T. A.,
        The MIT Press, Cambridge, Massachusetts, 2019.
    \endverbatim

SourceFiles
    gradientDescentMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef gradientDescentMinimizationScheme_H
#define gradientDescentMinimizationScheme_H

#include "minimizationScheme.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class gradientDescentMinimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class gradientDescentMinimizationScheme
:
    public minimizationScheme
{

    //- Name of the error
    virtual word errorName() const
    {
        return "delta";
    }


public:

    //- Runtime type information
    TypeName("gradientDescent");

    // Constructors

        //- Construct for given minimizationScheme
        gradientDescentMinimizationScheme
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );

    //- Destructor
    virtual ~gradientDescentMinimizationScheme()
    {}


    // Member Functions

        //- Solve the root
        virtual tmp<scalarField> minimize
        (
            const scalarField& x0,
            const scalarField& x1,
            const scalarField& x2,
            const label li
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

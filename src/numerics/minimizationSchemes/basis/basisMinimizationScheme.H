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
    Foam::minimizationSchemes::basis

Description
    Basic class for schemes using a basis line search

References:
    \verbatim
        "Algorithms for Optimization",
        Kochenderfer, M. J., Wheeler, T. A.,
        The MIT Press, Cambridge, Massachusetts, 2019.
    \endverbatim

SourceFiles
    basisMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef basisMinimizationScheme_H
#define basisMinimizationScheme_H

#include "minimizationScheme.H"
#include "cmptLineSearchEquation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace minimizationSchemes
{

/*---------------------------------------------------------------------------*\
                    Class basis Declaration
\*---------------------------------------------------------------------------*/

class basis
:
    public minimizationScheme
{
protected:
    // Protected data

        //- Component line search equation
        mutable cmptLineSearchEquation cmptLsEqn_;


    // Protected functions

        //- Search in the given direction
        void searchDir
        (
            const scalarList& x0,
            const label li,
            const label dir,
            const label sign,
            scalarList& xNew
        ) const;


public:

    //- Runtime type information
    TypeName("basis");

    // Constructors

        //- Construct for given minimizationScheme
        basis
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );

    //- Destructor
    virtual ~basis()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace minimizationSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

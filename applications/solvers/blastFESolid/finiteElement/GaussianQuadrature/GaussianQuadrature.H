/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
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
    Foam::GaussianQuadrature

Description
    Quadrature calculations for Legendre and Lobatto

SourceFiles
    GaussianQuadrature.C

\*---------------------------------------------------------------------------*/

#ifndef GaussianQuadrature_H
#define GaussianQuadrature_H

#include "scalar.H"
#include "List.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class GaussianQuadrature Declaration
\*---------------------------------------------------------------------------*/

class GaussianQuadrature
{
public:

    // Binomial coefficient
    static scalar binCoeff(const scalar n, const scalar k);

    //- Legendre polynomial
    static void calcLegendre
    (
        const label n,
        List<scalar>& x,
        List<scalar>& w
    );

    //- Lobatto polynomial
    static void calcLobatto
    (
        const label n,
        List<scalar>& x,
        List<scalar>& w
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // GaussianQuadrature_H

// ************************************************************************* //

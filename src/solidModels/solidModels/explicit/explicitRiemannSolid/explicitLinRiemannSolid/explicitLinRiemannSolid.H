/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    explicitLinRiemannSolid

Description
    Non-linear variation of the explicitTotalLag solid
    Modified version of the explicit solid dynamics library.

    Modifications were made to allow for the use of more complicated
    boundary conditions enforced through displacement, and additional
    generalization


Reference:

    Haider, J., Lee, C. H., Gill, A. J., Bonet, J. (2017),
    "A first-order hyperbolic framework for large strain computational
    solid dynamics: An upwind cell centred Total Lagrangian scheme",
    International Journal for Numerical Methods in Engineering,
    109, 407-456.

    https://github.com/jibranhaider/explicitSolidDynamics

SourceFiles
    explicitLinGeomTotalDispSolid.C

Also see:
    explicitLinGeomTotalDispSolid.C
    explicitTotalLagSolid.C

\*---------------------------------------------------------------------------*/

#ifndef explicitLinRiemannSolid_H
#define explicitLinRiemannSolid_H

#include "explicitRiemannSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class explicitLinRiemannSolid Declaration
\*---------------------------------------------------------------------------*/

class explicitLinRiemannSolid
:
    public explicitRiemannSolid
{
public:

    //- Runtime type information
    TypeName("explicitLinearRiemann");

    // Constructors

        //- Construct from components
        explicitLinRiemannSolid(dynamicFvMesh& mesh);

    // Destructor

        virtual ~explicitLinRiemannSolid()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

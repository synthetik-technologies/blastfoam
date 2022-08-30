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
    Foam::HRIC

Description
    HRIC differencing scheme

Reference
    Muzaferija, S. and Peric, M:
    Computation of free surface flows using interface-tracking and
    interface-capturing methods
    Chapter 2 in Mahrenholtz, O. and Markiewicz, M. (eds):
    Nonlinear Water Wate Interaction,
    Computational Mechanics Publications, Southampton, 1998.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    HRIC.C

\*---------------------------------------------------------------------------*/

#ifndef HRIC_H
#define HRIC_H

#include "limitedSurfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class HRIC Declaration
\*---------------------------------------------------------------------------*/

class HRIC
:
    public limitedSurfaceInterpolationScheme<scalar>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        HRIC(const HRIC&);

        //- Disallow default bitwise assignment
        void operator=(const HRIC&);


        //- Calculate individual face value
        scalar phiface
        (
            const scalar cdWeight,
            const scalar faceFlux,
            const scalar& phiP,
            const scalar& phiN,
            const vector& gradcP,
            const vector& gradcN,
            const vector d
        ) const;

        //- Calculate individual limiter
        scalar limiter
        (
            const scalar cdWeight,
            const scalar faceFlux,
            const scalar& phiP,
            const scalar& phiN,
            const vector& gradcP,
            const vector& gradcN,
            const scalar Cof,
            const vector d
        ) const;

        //- Calculate individual weight
        scalar weight
        (
            const scalar cdWeight,
            const scalar faceFlux,
            const scalar& phiP,
            const scalar& phiN,
            const vector& gradcP,
            const vector& gradcN,
            const scalar Cof,
            const vector d
        ) const;

public:

    //- Runtime type information
    TypeName("HRIC");


    // Constructors

        //- Construct from mesh and faceFlux
        HRIC
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            limitedSurfaceInterpolationScheme<scalar>(mesh, faceFlux)
        {}

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        HRIC
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<scalar>(mesh, is)
        {}

        //- Construct from mesh, faceFlux and Istream
        HRIC
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream&
        )
        :
            limitedSurfaceInterpolationScheme<scalar>(mesh, faceFlux)
        {}


    // Destructor

        virtual ~HRIC()
        {}


    // Member Functions

        //- Return the interpolation limiter
        virtual tmp<surfaceScalarField> limiter
        (
            const volScalarField&
        ) const;

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const volScalarField&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

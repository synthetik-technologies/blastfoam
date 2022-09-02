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
    Foam::mHRICLimiter

Description
    Class with limiter function which returns the limiter for the
    modified HRIC differencing scheme

    Used in conjunction with the template class LimitedScheme.

Author
    Heebum Lee, Seoul National University
    Rewrite and Co number sensitisation by Hrvoje Jasak, Wikki Ltd.

SourceFiles
    mHRIC.C

\*---------------------------------------------------------------------------*/

#ifndef mHRIC_H
#define mHRIC_H

#include "limitedSurfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mHRIC Declaration
\*---------------------------------------------------------------------------*/

class mHRIC
:
    public limitedSurfaceInterpolationScheme<scalar>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        mHRIC(const mHRIC&);

        //- Disallow default bitwise assignment
        void operator=(const mHRIC&);


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
    TypeName("mHRIC");


    // Constructors

        //- Construct from mesh and faceFlux
        mHRIC
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
        mHRIC
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<scalar>(mesh, is)
        {}

        //- Construct from mesh, faceFlux and Istream
        mHRIC
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream&
        )
        :
            limitedSurfaceInterpolationScheme<scalar>(mesh, faceFlux)
        {}


    // Destructor

        virtual ~mHRIC()
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

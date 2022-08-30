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
    Foam::CICSAM

Description
    CICSAM differencing scheme.

        Blending factor range 0 -> 1  STRICT !

        0   : Ultimate-Quickest (UQ)
        1   : Hyper-C           (HC)
        0.5 : recommended value for CICSAM

    Recommended changes to Blending factor :

        If interface becomes wavy use more UQ by reducing blending factor
        This may cause unacceptable interface smearing.

        If unconditional compression is required then use HC.  This
        may cause unacceptable wavyness.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    CICSAM.C

\*---------------------------------------------------------------------------*/

#ifndef CICSAM_H
#define CICSAM_H

#include "limitedSurfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class CICSAM Declaration
\*---------------------------------------------------------------------------*/

class CICSAM
:
    public limitedSurfaceInterpolationScheme<scalar>
{
    // Private data

        //- Scheme factor
        const scalar k_;


    // Static private data

        //- Lower bound
        static const scalar lowerBound_;

        //- Upper bound
        static const scalar upperBound_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CICSAM(const CICSAM&);

        //- Disallow default bitwise assignment
        void operator=(const CICSAM&);


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
    TypeName("CICSAM");


    // Constructors

        //- Construct from mesh, faceFlux and blendingFactor
        CICSAM
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const scalar k
        )
        :
            limitedSurfaceInterpolationScheme<scalar>(mesh, faceFlux),
            k_(k)
        {}

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        CICSAM
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<scalar>(mesh, is),
            k_(readScalar(is))
        {
            if (k_ < 0 || k_ > 1)
            {
                FatalIOErrorIn("CICSAMWeight(Istream& is)", is)
                    << "coefficient = " << k_
                    << " should be >= 0 and <= 1"
                    << exit(FatalIOError);
            }
        }

        //- Construct from mesh, faceFlux and Istream
        CICSAM
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<scalar>(mesh, faceFlux),
            k_(readScalar(is))
        {
            if (k_ < 0 || k_ > 1)
            {
                FatalIOErrorIn("CICSAMWeight(Istream& is)", is)
                    << "coefficient = " << k_
                    << " should be >= 0 and <= 1"
                    << exit(FatalIOError);
            }
        }


    // Destructor

        virtual ~CICSAM()
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

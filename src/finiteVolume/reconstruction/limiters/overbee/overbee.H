/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
     \\/     M anipulation  | Synthetik Applied Technology
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
    Foam::overbeeLimiter

Description
    Class with limiter function which returns the limiter for the
    overbee differencing scheme based on r obtained from the LimiterFunc
    class.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    overbee.C

\*---------------------------------------------------------------------------*/

#ifndef overbeeLimiter_H
#define overbeeLimiter_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class overbeeLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class overbeeLimiter
:
    public LimiterFunc
{
    const scalar beta_;

public:

    overbeeLimiter(Istream& is)
    :
        beta_(readScalar(is))
    {
        if (beta_ < 1 || beta_ > 2)
        {
            FatalIOErrorInFunction(is)
                << "coefficient = " << beta_
                << " should be >= 1 and <= 2"
                << exit(FatalIOError);
        }
    }

    scalar limiter
    (
        const scalar,
        const scalar faceFlux,
        const typename LimiterFunc::phiType& phiP,
        const typename LimiterFunc::phiType& phiN,
        const typename LimiterFunc::gradPhiType& gradcP,
        const typename LimiterFunc::gradPhiType& gradcN,
        const vector& d
    ) const
    {
        scalar r = LimiterFunc::r
        (
            faceFlux, phiP, phiN, gradcP, gradcN, d
        );

        return max
        (
            0.0,
            min
            (
                2.0,
                min
                (
                    2.0*r,
                    max
                    (
                        min(2.0*r, beta_),
                        min((2.0 - beta_)*r + 2.0*(beta_ - 1.0), r)
                    )
                )
            )
        );
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    BarthJespersenLimiter

Description
    Barth-Jespersen limiter

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef BarthJaspersenGradientLimiter_H
#define BarthJaspersenGradientLimiter_H

#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace gradientLimiters
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class BarthJaspersen
{

public:

    // Constructors

        BarthJaspersen(Istream& schemeData)
        {}


    // Member Functions

        //- Set scalar limiter value
        inline void limitFaceCmpt
        (
            scalar& lim,
            const scalar& cellVolume,
            const scalar& deltaOneMax,
            const scalar& deltaOneMin,
            const scalar& extrapolate
        ) const
        {
            if (mag(extrapolate) < SMALL)
            {
                return;
            }

            if (extrapolate - deltaOneMax > SMALL)
            {
                lim = min(lim, deltaOneMax/extrapolate);
            }
            else if (extrapolate - deltaOneMin < -SMALL)
            {
                lim = min(lim, deltaOneMin/extrapolate);
            }
        }

        //- Set Type limiter
        template<class Type>
        inline void limitFace
        (
            Type& lim,
            const scalar& cellVolume,
            const Type& deltaOneMax,
            const Type& deltaOneMin,
            const Type& extrapolate
        ) const
        {
            for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
            {
                limitFaceCmpt
                (
                    setComponent(lim, cmpt),
                    cellVolume,
                    component(deltaOneMax, cmpt),
                    component(deltaOneMin, cmpt),
                    component(extrapolate, cmpt)
                );
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace gradientLimiters
} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

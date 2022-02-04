/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    nonLinearGeometry

Description
    Enumerator inidacting whether the mathematical is employs a linear geometry,
    total Lagrangian nonlinear geometry or updated Lagrangian nonlinear geometry
    approach.

Author
    Hrvoje Jasak
    Philip Cardiff, UCD.

SourceFiles
    nonLinearGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinearGeometry_H
#define nonLinearGeometry_H

#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class nonLinearGeometry Declaration
\*---------------------------------------------------------------------------*/

class nonLinearGeometry
{
public:

        //- Non-linear solver options
        enum nonLinearType
        {
            LINEAR_GEOMETRY,
            UPDATED_LAGRANGIAN,
            TOTAL_LAGRANGIAN
        };

        static const NamedEnum<nonLinearType, 3> nonLinearNames_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

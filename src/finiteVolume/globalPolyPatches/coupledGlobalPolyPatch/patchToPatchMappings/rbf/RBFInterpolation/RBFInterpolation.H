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
    RBFInterpolation

Description
    Radial basis function interpolation class

Description
    Interpolation class which uses Radial Basis Functions to interpolate
    field from given data points to arbitrary set of points.

    The coefficient vectors, alpha and beta are determined by solving
    the system:

    | db | = | Mbb Pb | | alpha |
    | 0  |   | Pb  0  | |  beta |

    where db are the given field values at data carrier points.
    Mbb the boundary RBF correlation matrix (NbxNb), containing RBF evaluations
    at the boundary nodes, and Pb some linear polynomial matrix (Nbx4).

    In cases where far field data is not of interest, a cutoff function
    is used to eliminate unnecessary data points in the far field

Author
    Frank Bos, TU Delft.  All rights reserved.
    Dubravko Matijasevic, FSB Zagreb.
    Reorganisation by Hrvoje Jasak, Wikki Ltd.

SourceFiles
    RBFInterpolation.C
    RBFInterpolationTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef RBFInterpolation_H
#define RBFInterpolation_H

#include "dictionary.H"
#include "RBFFunction.H"
#include "point.H"
#include "Switch.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class RBFInterpolation Declaration
\*---------------------------------------------------------------------------*/

class RBFInterpolation
{
    // Private data

        //- Reference to control points
        const vectorField& controlPoints_;

        //- Reference to all points
        const vectorField& dataPoints_;

        //- RBF function
        autoPtr<RBFFunction> RBF_;

        //- Interpolation matrix
        mutable scalarSquareMatrix* BPtr_;

        //- Focal point for cut-off radii
        point focalPoint_;

        //- Inner cut-off radius
        scalar innerRadius_;

        //- Outer cut-off radius
        scalar outerRadius_;

        //- Add polynomials to RBF matrix
        Switch polynomials_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const RBFInterpolation&);


        //- Return interpolation matrix
        const scalarSquareMatrix& B() const;

        //- Calculate interpolation matrix
        void calcB() const;

        //- Clear out
        void clearOut();


public:

    // Constructors

        //- Construct from components
        RBFInterpolation
        (
            const dictionary& dict,
            const vectorField& controlPoints,
            const vectorField& dataPoints
        );

        //- Construct from components
        RBFInterpolation
        (
            const word& type,
            const dictionary& dict,
            const vectorField& controlPoints,
            const vectorField& dataPoints
        );

        //- Construct as copy
        RBFInterpolation(const RBFInterpolation&);


    // Destructor

        ~RBFInterpolation();


    // Member Functions

        //- Return reference to control points
        const vectorField& controlPoints() const
        {
            return controlPoints_;
        }

        //- Reference to all points
        const vectorField& dataPoints() const
        {
            return dataPoints_;
        }

        //- Return focalPoint
        const point& focalPoint() const
        {
            return focalPoint_;
        }

        //- Access focalPoint
        point& focalPoint()
        {
            return focalPoint_;
        }

        //- Return inner radius
        scalar innerRadius() const
        {
            return innerRadius_;
        }

        //- Access inner radius
        scalar& innerRadius()
        {
            return innerRadius_;
        }

        //- Return outer radius
        scalar outerRadius() const
        {
            return outerRadius_;
        }

        //- Access outer radius
        scalar& outerRadius()
        {
            return outerRadius_;
        }

        //- Interpolate
        template<class Type>
        tmp<Field<Type>> interpolate(const Field<Type>& ctrlField) const;

        //- Interpolate
        template<class Type>
        void interpolate
        (
            const Field<Type>& ctrlField,
            Field<Type>& resField
        ) const;

        //- Move points
        void movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "RBFInterpolationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

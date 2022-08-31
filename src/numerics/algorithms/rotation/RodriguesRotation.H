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
    RodriguesRotation

Description
    Compute rotation tensor base on Rodrigues' Rotation Formula
    http://mathworld.wolfram.com/RodriguesRotationFormula.html
    http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula

Author
    Martin Beaudoin, Hydro-Quebec, (2008, 2011)

SourceFiles
    RodriguesRotation.C

\*---------------------------------------------------------------------------*/

#ifndef RodriguesRotation_H
#define RodriguesRotation_H

#include "tmp.H"
#include "vector.H"
#include "tensor.H"
#include "tensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Compute rotation tensor using Rodrigues' Rotation Formula
// Formulation taken from:
// http://mathworld.wolfram.com/RodriguesRotationFormula.html
//
// Requires: rotation axis + rotation angle
tensor RodriguesRotation
(
    const vector& rotationAxis,
    const scalar& rotationAngle,
    const bool inDegrees = true
);

// Compute rotation tensor using Rodrigues' Rotation Formula
// Slightly different formulation taken from:
// http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
//
// Requires: rotation axis + rotation angle
tensor RodriguesRotation2
(
    const vector& rotationAxis,
    const scalar& rotationAngle,
    const bool inDegrees = true
);

// Compute rotation tensor using Rodrigues' Rotation Formula
// http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
//
// Requires: rotation axis + initial vector + final vector

tensor RodriguesRotation
(
    const vector& rotationAxis,
    const vector& vi,
    const vector& vf
);

// Compute rotation tensorField using Rodrigues' Rotation Formula
// http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
//
// Requires: rotation axis + initial vectorField + final vectorField
tmp<tensorField> RodriguesRotation
(
    const vector& rotationAxis,
    const vectorField& vi,
    const vectorField& vf
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

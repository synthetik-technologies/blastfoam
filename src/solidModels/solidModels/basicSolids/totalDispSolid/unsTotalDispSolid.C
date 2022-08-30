/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "unsTotalDispSolid.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void unsTotalDispSolid::updateDisplacement()
{
    // Update the total displacement
    DD() = D() - D().oldTime();

    // Interpolate D to pointD
    mechanical().interpolate(D(), pointD(), false);

    // Increment of displacement
    pointDD() = pointD() - pointD().oldTime();

    // Update gradient of displacement
    mechanical().grad(D(), pointD(), gradD(), gradDf_);
    mechanical().grad(D(), pointD(), gradD());

    // Update gradient of total displacement
    gradDD() = gradD() - gradD().oldTime();
    gradDDf_ = gradDf_ - gradDf_.oldTime();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unsTotalDispSolid::unsTotalDispSolid
(
    const word& type,
    dynamicFvMesh& mesh,
    const nonLinearGeometry::nonLinearType nonLinear,
    const bool isSolid
)
:
    totalDispSolid(type, mesh, nonLinear, isSolid),
    sigmaf_
    (
        IOobject
        (
            "sigmaf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    gradDf_
    (
        IOobject
        (
            "grad(" + D().name() + ")f",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    gradDDf_
    (
        IOobject
        (
            "grad(" + DD().name() + ")f",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    )
{
    // For consistent restarts, we will calculate the gradient field
    mechanical().grad(D(), pointD(), gradD(), gradDf_);
    mechanical().grad(DD(), pointDD(), gradDD(), gradDDf_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

} // End namespace Foam

// ************************************************************************* //

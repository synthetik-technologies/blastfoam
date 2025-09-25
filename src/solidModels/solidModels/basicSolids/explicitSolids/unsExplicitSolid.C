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

\*---------------------------------------------------------------------------*/

#include "unsExplicitSolid.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * * //

void Foam::solidModels::unsExplicitSolid::updateDisplacement()
{
    // Update gradient of displacement increment
    this->mechanical().grad(this->D(), this->gradD());

    // Interpolate cell displacements to vertices
    this->mechanical().interpolate(this->DD(), this->pointDD());

    // Update gradient of displacement
    this->mechanical().grad
    (
        this->D(),
        this->pointD(),
        this->gradD(),
        this->gradDf()
    );

    // Update gradient of total displacement
    this->gradDD() = this->gradD() - this->gradD().oldTime();
    this->gradDDf() = this->gradDf() - this->gradDf().oldTime();

    // Increment of displacement
    this->pointD() == this->pointD().oldTime() - this->pointDD();
}


void Foam::solidModels::unsExplicitSolid::updateWavespeeds()
{
    wavespeed_ =
        fvc::interpolate(sqrt(this->mechanical().elasticModulus()/this->rho()));
    sWavespeed_ =
        fvc::interpolate
        (
            sqrt(this->mechanical().shearModulus()/this->rho())
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidModels::unsExplicitSolid::unsExplicitSolid
(
    const word& type,
    fvMesh& mesh,
    const nonLinearGeometry::nonLinearType nonLinear,
    const bool isSolid
)
:
    ExplicitSolidBase<UnsIncrementalSolid<unsDSolidModel>>
    (
        type,
        mesh,
        nonLinear,
        isSolid
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //

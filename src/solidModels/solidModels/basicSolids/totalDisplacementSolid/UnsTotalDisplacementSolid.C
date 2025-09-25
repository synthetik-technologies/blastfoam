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

#include "UnsTotalDisplacementSolid.H"

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

template<class UnsSolidModel>
void Foam::solidModels::UnsTotalDisplacementSolid<UnsSolidModel>::
updateDisplacement()
{
    // Update the total displacement
    this->DD() = this->D() - this->D().oldTime();

    // Interpolate D to pointD
    this->mechanical().interpolate(this->D(), this->pointD(), false);

    // Increment of displacement
    this->pointDD() = this->pointD() - this->pointD().oldTime();

    // Update gradient of displacement
    // this->mechanical().grad(this->D(), this->pointD(), this->gradD());
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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class UnsSolidModel>
Foam::solidModels::UnsTotalDisplacementSolid<UnsSolidModel>::
UnsTotalDisplacementSolid
(
    const word& type,
    fvMesh& mesh,
    const nonLinearGeometry::nonLinearType nonLinear,
    const bool isSolid
)
:
    UnsSolidModel(type, mesh, nonLinear, incremental(), isSolid),
    impKf_("impKf", this->mechanical().impKf())
{
    this->isRequired(this->D(), type);

    // Interpolate D to pointD
    this->mechanical().interpolate(this->D(), this->pointD(), false);

    // For consistent restarts, we will calculate the gradient field
    this->mechanical().grad
    (
        this->D(),
        this->pointD(),
        this->gradD(),
        this->gradDf()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //

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

#include "IncrementalSolid.H"

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

template<class SolidModel>
void Foam::solidModels::IncrementalSolid<SolidModel>::updateDisplacement()
{
    // Update the total displacement
    this->D() = this->D().oldTime() + this->DD();

    // Interpolate DD to pointDD
    this->mechanical().interpolate(this->DD(), this->pointDD(), false);

    // Update gradient of displacement increment
    this->mechanical().grad(this->DD(), this->gradDD());

    // Update gradient of total displacement
    this->gradD() = this->gradD().oldTime() + this->gradDD();

    // Total displacement at points
    this->pointD() == this->pointD().oldTime() + this->pointDD();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SolidModel>
Foam::solidModels::IncrementalSolid<SolidModel>::IncrementalSolid
(
    const word& type,
    fvMesh& mesh,
    const nonLinearGeometry::nonLinearType nonLinear,
    const bool isSolid
)
:
    SolidModel(type, mesh, nonLinear, incremental(), isSolid),
    impK_("impK", this->mechanical().impK()),
    impKf_("impKf", this->mechanical().impKf())
{
    this->isRequired(this->DD(), type);

    // For consistent restarts, we will calculate the gradient field
    this->mechanical().grad(this->D(), this->gradD());
    this->gradD().storeOldTimes();

    this->mechanical().grad(this->DD(), this->gradDD());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //

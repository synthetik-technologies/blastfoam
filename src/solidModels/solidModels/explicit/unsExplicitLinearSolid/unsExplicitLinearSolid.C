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

#include "unsExplicitLinearSolid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(unsExplicitLinearSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, unsExplicitLinearSolid, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unsExplicitLinearSolid::unsExplicitLinearSolid
(
    dynamicFvMesh& mesh
)
:
    UnsLinearGeomSolid<unsExplicitSolid>(typeName, mesh)
{
    a_.oldTime();
    U().oldTime();

    // Update stress
    update();

    // Update initial acceleration
    a_ = fvc::div(this->tractionSf())/rho();
    a_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void unsExplicitLinearSolid::update(const bool correctSigma)
{
    UnsLinearGeomSolid<unsExplicitSolid>::update(correctSigma);

    if (correctSigma)
    {
        this->impKf_ = this->mechanical().impKf();
        wavespeed_ = sqrt(this->impKf_/fvc::interpolate(rho()));
    }
}


tmp<surfaceVectorField> unsExplicitLinearSolid::tractionSf() const
{
    return this->mesh().Sf() & this->sigmaf();
}


tmp<volVectorField> unsExplicitLinearSolid::divStress() const
{
    return fvc::div(this->mesh().Sf() & this->sigmaf());
}


bool unsExplicitLinearSolid::evolve()
{
    this->solveMomentum();
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //

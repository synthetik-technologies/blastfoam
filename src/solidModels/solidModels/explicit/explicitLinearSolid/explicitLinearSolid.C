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

#include "explicitLinearSolid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitLinearSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, explicitLinearSolid, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitLinearSolid::explicitLinearSolid(fvMesh& mesh)
:
    LinearGeomSolid<explicitSolid>(typeName, mesh)
{
    U().oldTime();

    // Update stress
    update();

    // Update initial acceleration
    a_ = fvc::div(sigma(), "div(sigma)")/rho();
    a_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void explicitLinearSolid::update(const bool correctSigma)
{
    // Update displacements and gradients
    LinearGeomSolid<explicitSolid>::update(correctSigma);

    if (correctSigma)
    {
        this->impKf_ = this->mechanical().impKf();
        this->updateWavespeeds();
        this->wavespeed_ = sqrt(this->impKf_/fvc::interpolate(this->rho()));
    }
}


tmp<surfaceVectorField> explicitLinearSolid::tractionSf() const
{
    return fvc::dotInterpolate(this->mesh().Sf(), this->sigma());
}


tmp<volVectorField> explicitLinearSolid::divStress() const
{
    return fvc::div(this->sigma(), "div(sigma)");
}


bool explicitLinearSolid::evolve()
{
    this->solveMomentum();
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //

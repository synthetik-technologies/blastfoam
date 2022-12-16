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

#include "explicitTotalLagrangianSolid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitTotalLagrangianSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, explicitTotalLagrangianSolid, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitTotalLagrangianSolid::explicitTotalLagrangianSolid
(
    dynamicFvMesh& mesh
)
:
    TotalLagrangianGeomSolid<explicitNonLinearSolid>(typeName, mesh)
{
    // Update stress
    update();

    // Update initial acceleration
    a_ = fvc::div(sigma(), "div(sigma)")/rho();
    a_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void explicitTotalLagrangianSolid::update(const bool correctSigma)
{
    TotalLagrangianGeomSolid<explicitNonLinearSolid>::update(correctSigma);

    if (correctSigma)
    {
        this->impKf_ = this->mechanical().impKf();
        wavespeed_ = sqrt(this->impKf_/fvc::interpolate(rho()));
    }
    this->updateWavespeeds();
}


tmp<volVectorField> explicitTotalLagrangianSolid::divStress() const
{
    return fvc::div(this->P(), "div(sigma)");
}


tmp<surfaceVectorField> explicitTotalLagrangianSolid::tractionSf() const
{
    return fvc::dotInterpolate(this->mesh().Sf(), this->P());
}


bool explicitTotalLagrangianSolid::evolve()
{
    this->solveMomentum();
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //

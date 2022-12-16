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

#include "explicitUpdatedLagrangianSolid.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitUpdatedLagrangianSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, explicitUpdatedLagrangianSolid, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitUpdatedLagrangianSolid::explicitUpdatedLagrangianSolid
(
    dynamicFvMesh& mesh
)
:
    UpdatedLagrangianGeomSolid<explicitNonLinearSolid>(typeName, mesh)
{
    a_.oldTime();
    U().oldTime();

    // Update stress
    update();

    // Update initial acceleration
    if (!a_.headerOk())
    {
        a_ = fvc::div(sigma(), "div(sigma)")/rho();
        a_.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void explicitUpdatedLagrangianSolid::update(const bool correctSigma)
{
    UpdatedLagrangianGeomSolid<explicitNonLinearSolid>::update(correctSigma);

    if (correctSigma)
    {
        updateWavespeeds();
        this->updateWavespeeds();
    }
}


tmp<surfaceVectorField> explicitUpdatedLagrangianSolid::tractionSf() const
{
    return fvc::dotInterpolate(this->mesh().Sf(), this->P());
}


tmp<volVectorField> explicitUpdatedLagrangianSolid::divStress() const
{
    return fvc::div(this->P(), "div(sigma)");
}


bool explicitUpdatedLagrangianSolid::evolve()
{
    this->solveMomentum();
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //

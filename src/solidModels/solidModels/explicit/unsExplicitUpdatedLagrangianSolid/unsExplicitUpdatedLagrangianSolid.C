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

#include "unsExplicitUpdatedLagrangianSolid.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(unsExplicitUpdatedLagrangianSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, unsExplicitUpdatedLagrangianSolid, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unsExplicitUpdatedLagrangianSolid::unsExplicitUpdatedLagrangianSolid
(
    fvMesh& mesh
)
:
    UnsUpdatedLagrangianGeomSolid<unsExplicitNonLinearSolid>(typeName, mesh)
{
    a_.oldTime();
    U().oldTime();

    // Update stress
    update();

    // Update initial acceleration
    a_.primitiveFieldRef() =
        fvc::div(sigma(), "div(sigma)")().internalField()
       /(rho().internalField());
    a_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void unsExplicitUpdatedLagrangianSolid::update(const bool correctSigma)
{
    UnsUpdatedLagrangianGeomSolid<unsExplicitNonLinearSolid>::update(correctSigma);

    if (correctSigma)
    {
        this->sigma() = fvc::average(this->sigmaf());
        this->updateWavespeeds();
        this->impKf_ = this->mechanical().impKf();
    }
}


tmp<surfaceVectorField> unsExplicitUpdatedLagrangianSolid::tractionSf() const
{
    return this->mesh().Sf() & this->Pf();
}


tmp<volVectorField> unsExplicitUpdatedLagrangianSolid::divStress() const
{
    return fvc::div(this->mesh().Sf() & this->Pf());
}


bool unsExplicitUpdatedLagrangianSolid::evolve()
{
    this->solveMomentum();
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //

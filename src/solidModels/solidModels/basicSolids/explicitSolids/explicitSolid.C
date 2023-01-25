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

#include "explicitSolid.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * * //

void Foam::solidModels::explicitSolid::updateDisplacement()
{
    // Update gradient of displacement increment
    this->mechanical().grad(this->D(), this->gradD());
    this->mechanical().grad(this->DD(), this->gradDD());

    // Update gradient of total displacement
    this->gradDD() = this->gradD() - this->gradD().oldTime();

    // Interpolate cell displacements to vertices
    this->mechanical().interpolate(this->D(), this->pointD());

    // Increment of displacement
    this->pointDD() = this->pointD() - this->pointD().oldTime();
}


void Foam::solidModels::explicitSolid::correctUBCs
(
    volVectorField& U
)
{
    // U.correctBoundaryConditions();
    // volVectorField::Boundary& bU = U.boundaryFieldRef();
    // const volVectorField::Boundary& bD = this->solutionD().boundaryField();
    // forAll(bU, patchi)
    // {
    //     if (isA<solidTractionFvPatchVectorField>(bD[patchi]))
    //     {
    //         const fvPatch& patch = this->mesh().boundary()[patchi];
    //
    //         fvPatchVectorField& pU = bU[patchi];
    //         const solidTractionFvPatchVectorField& pD =
    //             dynamicCast<const solidTractionFvPatchVectorField>(bD[patchi]);
    //         vectorField n(this->nf(patch));
    //         symmTensorField psigma(this->sigma(patch));
    //         tensorField nn(n*n);
    //
    //     tensorField St
    //     (
    //         nn/wavespeed_.boundaryField()[patchi]
    //       + (I - nn)/sWavespeed_.boundaryField()[patchi]
    //     );
    //
    //     pU =
    //         pU.internalField()
    //       + (
    //             St & (pD.traction() - pD.pressure()*n - (n & psigma))
    //         )/this->rho().boundaryField()[patchi];
    //     }
    // }
}


void Foam::solidModels::explicitSolid::updateWavespeeds()
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

Foam::solidModels::explicitSolid::explicitSolid
(
    const word& type,
    dynamicFvMesh& mesh,
    const nonLinearGeometry::nonLinearType nonLinear,
    const bool isSolid
)
:
    ExplicitSolidBase<IncrementalSolid<solidModel>>
    (
        type,
        mesh,
        nonLinear,
        isSolid
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //

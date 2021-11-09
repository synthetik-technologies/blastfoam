/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "solidTractionCoupledImmersedFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidTractionCoupledImmersedFvPatchVectorField::
solidTractionCoupledImmersedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    name_
    (
        p.boundaryMesh().mesh().dbDir()
    ),
    mapper_
    (
        p.boundaryMesh().mesh().time().lookupObject<immersedMeshMapper>
        (
            IOobject::groupName("immersedMeshMapper", name_)
        )
    ),
    ibm_(mapper_.immersedObject()),
    pRef_(0.0)
{}


solidTractionCoupledImmersedFvPatchVectorField::
solidTractionCoupledImmersedFvPatchVectorField
(
    const solidTractionCoupledImmersedFvPatchVectorField& psf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(psf, p, iF, mapper),
    name_
    (
        p.boundaryMesh().mesh().dbDir()
    ),
    mapper_
    (
        p.boundaryMesh().mesh().time().lookupObject<immersedMeshMapper>
        (
            IOobject::groupName("immersedMeshMapper", name_)
        )
    ),
    ibm_(mapper_.immersedObject()),
    pRef_(psf.pRef_)
{}


solidTractionCoupledImmersedFvPatchVectorField::
solidTractionCoupledImmersedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF, dict),
    name_
    (
        p.boundaryMesh().mesh().dbDir()
    ),
    mapper_
    (
        p.boundaryMesh().mesh().time().lookupObject<immersedMeshMapper>
        (
            IOobject::groupName("immersedMeshMapper", name_)
        )
    ),
    ibm_(mapper_.immersedObject()),
    pRef_(dict.lookupOrDefault<scalar>("pRef", 0.0))
{}


solidTractionCoupledImmersedFvPatchVectorField::
solidTractionCoupledImmersedFvPatchVectorField
(
    const solidTractionCoupledImmersedFvPatchVectorField& psf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(psf, iF),
    name_(psf.name_),
    mapper_(psf.mapper_),
    ibm_(psf.ibm_),
    pRef_(psf.pRef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidTractionCoupledImmersedFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (ibm_.pMesh().foundObject<volScalarField>("p"))
    {
        const volScalarField& volNbrp =
            ibm_.pMesh().lookupObject<volScalarField>("p");
        this->pressure() =
            mapper_.mapObjectToBoundary
            (
                ibm_.interpolateTo(volNbrp)()
            ) - pRef_;
    }
    else
    {
        this->pressure() = Zero;
    }

    if (ibm_.pMesh().foundObject<volSymmTensorField>("devTau"))
    {
        const volSymmTensorField& volNbrDevTau =
            ibm_.pMesh().lookupObject<volSymmTensorField>("devTau");
        this->traction() =
            mapper_.mapObjectToBoundary
            (
                ibm_.interpolateTo(volNbrDevTau)()
            ) & patch().Sf();
    }
    else
    {
        this->traction() = Zero;
    }

    solidTractionFvPatchVectorField::updateCoeffs();
}


void solidTractionCoupledImmersedFvPatchVectorField::write
(
    Ostream& os
) const
{
    solidTractionFvPatchVectorField::write(os);
    writeEntry(os, "pRef", pRef_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    solidTractionCoupledImmersedFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //

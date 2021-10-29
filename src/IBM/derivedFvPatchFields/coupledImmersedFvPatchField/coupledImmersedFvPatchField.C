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

#include "coupledImmersedFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
coupledImmersedFvPatchField<Type>::coupledImmersedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
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
    ibm_(mapper_.immersedObject())
{}


template<class Type>
coupledImmersedFvPatchField<Type>::coupledImmersedFvPatchField
(
    const coupledImmersedFvPatchField& psf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(p, iF),
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
    ibm_(mapper_.immersedObject())
{}


template<class Type>
coupledImmersedFvPatchField<Type>::coupledImmersedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
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
    ibm_(mapper_.immersedObject())
{}


template<class Type>
coupledImmersedFvPatchField<Type>::coupledImmersedFvPatchField
(
    const coupledImmersedFvPatchField& psf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(psf, iF),
    name_(psf.name_),
    mapper_(psf.mapper_),
    ibm_(psf.ibm_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void coupledImmersedFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    const FieldType& volNbr =
        ibm_.pMesh().lookupObject<FieldType>(this->internalField().name());
    this->operator==
    (
        mapper_.mapObjectToBoundary
        (
            ibm_.interpolateTo(volNbr)()
        )
    );

    fixedValueFvPatchField<Type>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //

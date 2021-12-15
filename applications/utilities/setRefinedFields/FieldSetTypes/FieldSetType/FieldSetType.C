/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "FieldSetType.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::FieldSetType<Type, Patch, Mesh>::FieldSetType
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& fieldName,
    const labelList& selectedIndices,
    Istream& is,
    const bool write
)
:
    mesh_(mesh),
    dict_(dict),
    fieldPtr_(lookupOrRead(fieldName)),
    selectedIndices_(selectedIndices),
    write_(write),
    good_(fieldPtr_.valid())
{}


template<class Type, template<class> class Patch, class Mesh>
Foam::FieldSetType<Type, Patch, Mesh>::FieldSetType
(
    const fvMesh& mesh,
    const dictionary& dict,
    const labelList& selectedIndices
)
:
    mesh_(mesh),
    dict_(dict),
    fieldPtr_(nullptr),
    selectedIndices_(selectedIndices),
    write_(false),
    good_(false)
{}


template<class Type>
Foam::VolFieldSetType<Type>::VolFieldSetType
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& fieldName,
    const labelList& selectedIndices,
    Istream& is,
    const bool write
)
:
    FieldSetType<Type, fvPatchField, volMesh>
    (
        mesh,
        dict,
        fieldName,
        selectedIndices,
        is,
        write
    )
{}


template<class Type>
Foam::VolFieldSetType<Type>::VolFieldSetType
(
    const fvMesh& mesh,
    const dictionary& dict,
    const labelList& selectedIndices
)
:
    FieldSetType<Type, fvPatchField, volMesh>(mesh, dict, selectedIndices)
{}


template<class Type>
Foam::SurfaceFieldSetType<Type>::SurfaceFieldSetType
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& fieldName,
    const labelList& selectedIndices,
    Istream& is,
    const bool write
)
:
    FieldSetType<Type, fvsPatchField, surfaceMesh>
    (
        mesh,
        dict,
        fieldName,
        selectedIndices,
        is,
        write
    )
{}


template<class Type>
Foam::SurfaceFieldSetType<Type>::SurfaceFieldSetType
(
    const fvMesh& mesh,
    const dictionary& dict,
    const labelList& selectedIndices
)
:
    FieldSetType<Type, fvsPatchField, surfaceMesh>(mesh, dict, selectedIndices)
{}


template<class Type>
Foam::PointFieldSetType<Type>::PointFieldSetType
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& fieldName,
    const labelList& selectedIndices,
    Istream& is,
    const bool write
)
:
    FieldSetType<Type, pointPatchField, pointMesh>
    (
        mesh,
        dict,
        fieldName,
        selectedIndices,
        is,
        write
    )
{}


template<class Type>
Foam::PointFieldSetType<Type>::PointFieldSetType
(
    const fvMesh& mesh,
    const dictionary& dict,
    const labelList& selectedIndices
)
:
    FieldSetType<Type, pointPatchField, pointMesh>(mesh, dict, selectedIndices)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::FieldSetType<Type, Patch, Mesh>::~FieldSetType()
{}


template<class Type>
Foam::VolFieldSetType<Type>::~VolFieldSetType()
{}


template<class Type>
Foam::SurfaceFieldSetType<Type>::~SurfaceFieldSetType()
{}


template<class Type>
Foam::PointFieldSetType<Type>::~PointFieldSetType()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
const Foam::fvMesh& Foam::FieldSetType<Type, Patch, Mesh>::mesh
(
    const UautoPtr<GeometricField<Type, fvPatchField, volMesh>>&
) const
{
    return mesh_;
}


template<class Type, template<class> class Patch, class Mesh>
const Foam::fvMesh& Foam::FieldSetType<Type, Patch, Mesh>::mesh
(
    const UautoPtr<GeometricField<Type, fvsPatchField, surfaceMesh>>&
) const
{
    return mesh_;
}


template<class Type, template<class> class Patch, class Mesh>
const Foam::pointMesh& Foam::FieldSetType<Type, Patch, Mesh>::mesh
(
    const UautoPtr<GeometricField<Type, pointPatchField, pointMesh>>&
) const
{
    return pointMesh::New(mesh_);
}


template<class Type, template<class> class Patch, class Mesh>
typename Foam::FieldSetType<Type, Patch, Mesh>::FieldType*
Foam::FieldSetType<Type, Patch, Mesh>::lookupOrRead(const word& fieldName) const
{
    if (mesh_.foundObject<FieldType>(fieldName))
    {
        return &mesh_.lookupObjectRef<FieldType>(fieldName);
    }


    // Check the current time directory
    IOobject fieldHeader
    (
        fieldName,
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ
    );

    // Check the "constant" directory
    if (!fieldHeader.typeHeaderOk<FieldType>(true))
    {
        fieldHeader = IOobject
        (
            fieldName,
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ
        );
    }

    // Check field exists
    if (fieldHeader.typeHeaderOk<FieldType>(true))
    {
        FieldType* fPtr
        (
            new FieldType(fieldHeader, mesh(fieldPtr_))
        );
        fPtr->store(fPtr);
        return &mesh_.lookupObjectRef<FieldType>(fieldName);
    }
    return nullptr;
}


template<class Type, template<class> class Patch, class Mesh>
void Foam::FieldSetType<Type, Patch, Mesh>::getInternalField
(
    const labelList& indices,
    const UIndirectList<vector>& pts,
    UIndirectList<Type>& f
)
{}


template<class Type, template<class> class Patch, class Mesh>
void Foam::FieldSetType<Type, Patch, Mesh>::getBoundaryField
(
    const label patchi,
    const labelList& indices,
    const UIndirectList<vector>& pts,
    UIndirectList<Type>& f
)
{}


template<class Type>
void Foam::VolFieldSetType<Type>::setField()
{
    const UIndirectList<vector> CInt(this->mesh_.C(), this->selectedIndices_);
    UIndirectList<Type> fInt(this->fieldPtr_(), this->selectedIndices_);
    this->getInternalField(this->selectedIndices_, CInt, fInt);

    typename GeometricField<Type, fvPatchField, volMesh>::
        Boundary& fieldBf = this->fieldPtr_->boundaryFieldRef();
    forAll(this->fieldPtr_->boundaryField(), patchi)
    {
        fieldBf[patchi] = fieldBf[patchi].patchInternalField();
    }

    if (this->write_)
    {
        if (!this->fieldPtr_->write())
        {
            FatalErrorInFunction
                << "Failed writing field " << this->fieldPtr_->name() << endl
                << exit(FatalError);
        }
    }
}

template<class Type>
void Foam::SurfaceFieldSetType<Type>::setField()
{
    const UIndirectList<vector> CfInt(this->mesh_.Cf(), this->selectedIndices_);
    UIndirectList<Type> fInt(this->fieldPtr_(), this->selectedIndices_);
    this->getInternalField(this->selectedIndices_, CfInt, fInt);

    if (this->write_)
    {
        if (!this->fieldPtr_->write())
        {
            FatalErrorInFunction
                << "Failed writing field " << this->fieldPtr_->name() << endl
                << exit(FatalError);
        }
    }
}

template<class Type>
void Foam::PointFieldSetType<Type>::setField()
{
    const UIndirectList<vector> pts(this->mesh_.points(), this->selectedIndices_);
    UIndirectList<Type> fInt(this->fieldPtr_(), this->selectedIndices_);
    this->getInternalField(this->selectedIndices_, pts, fInt);

    if (this->write_)
    {
        if (!this->fieldPtr_->write())
        {
            FatalErrorInFunction
                << "Failed writing field " << this->fieldPtr_->name() << endl
                << exit(FatalError);
        }
    }
}

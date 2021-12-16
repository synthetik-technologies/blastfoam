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
#include "valuePointPatchField.H"
#include "fieldSetOptions.H"

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
    noInternal_(false),
    write_(write),
    good_(fieldPtr_.valid())
{
    token t;
    while (true)
    {
        t = token(is);
        if (!t.isPunctuation())
        {
            break;
        }
        if (t.pToken() != '-')
        {
            break;
        }
        fieldSetOptions::options opt = fieldSetOptions::optionNames[word(is)];
        switch (opt)
        {
            case fieldSetOptions::SetBoundaries:
                boundaries_ = wordHashSet(is);
                break;
            case fieldSetOptions::SetAllBoundaries:
                boundaries_ = wordHashSet(mesh_.boundaryMesh().names());
                break;
            case fieldSetOptions::NoInternal:
                noInternal_ = true;
                break;
            default:
                FatalErrorInFunction
                    << "Unknown option \"" << opt << "\""<< nl
                    << "valid options are " << nl
                    << fieldSetOptions::optionNames.toc() << endl
                    << abort(FatalError);
        }
        if (noInternal_ && !boundaries_.size())
        {
            FatalErrorInFunction
                << "Internal field is not set and no boundaries were "
                << "specified" << endl
                << abort(FatalError);
        }
    }
    is.putBack(t);

}


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
Foam::tmp<Foam::Field<Type>> Foam::VolFieldSetType<Type>::getBoundary
(
    const label patchi,
    const GeometricField<Type, fvPatchField, volMesh>& f
) const
{
    return f.boundaryField()[patchi];
}


template<class Type>
void Foam::VolFieldSetType<Type>::setField()
{
    if (!this->noInternal_)
    {
        const UIndirectList<vector> CInt(this->mesh_.C(), this->selectedIndices_);
        UIndirectList<Type> fInt(this->fieldPtr_(), this->selectedIndices_);
        this->getInternalField(this->selectedIndices_, CInt, fInt);
    }

    labelHashSet cells(this->selectedIndices_);

    typename GeometricField<Type, fvPatchField, volMesh>::
        Boundary& fieldBf = this->fieldPtr_->boundaryFieldRef();
    forAll(this->fieldPtr_->boundaryField(), patchi)
    {
        const polyPatch& p = this->mesh_.boundaryMesh()[patchi];
        if (this->boundaries_.found(p.name()) && fieldBf[patchi].fixesValue())
        {
            const labelList& fCells = p.faceCells();
            labelHashSet faces;
            forAll(fCells, fi)
            {
                if (cells.found(fCells[fi]))
                {
                    faces.insert(fi);
                }
            }
            labelList indices(faces.toc());
            const UIndirectList<vector> pC
            (
                this->mesh_.C().boundaryField()[patchi],
                indices
            );
            UIndirectList<Type> pf(fieldBf[patchi], indices);
            this->getBoundaryField(patchi, indices, pC, pf);
        }
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
Foam::tmp<Foam::Field<Type>> Foam::SurfaceFieldSetType<Type>::getBoundary
(
    const label patchi,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& f
) const
{
    return f.boundaryField()[patchi];
}


template<class Type>
void Foam::SurfaceFieldSetType<Type>::setField()
{
    labelHashSet faces(this->selectedIndices_);
    labelList indices(this->selectedIndices_);
    label I = 0;
    forAll(indices, i)
    {
        if (indices[i] < this->mesh_.nInternalFaces())
        {
            indices[I++] = indices[i];
        }
    }
    indices.resize(I);

    if (!this->noInternal_)
    {
        const UIndirectList<vector> CfInt(this->mesh_.Cf(), indices);
        UIndirectList<Type> fInt(this->fieldPtr_(), indices);
        this->getInternalField(indices, CfInt, fInt);
    }

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
        Boundary& fieldBf = this->fieldPtr_->boundaryFieldRef();

    forAll(this->fieldPtr_->boundaryField(), patchi)
    {
        const polyPatch& p = this->mesh_.boundaryMesh()[patchi];
        if (this->boundaries_.found(p.name()) && fieldBf[patchi].fixesValue())
        {
            indices = identity(p.size());
            I = 0;
            forAll(p, fi)
            {
                if (faces.found(fi + p.start()))
                {
                    indices[I++] = fi;
                }
            }
            if (!I)
            {
                continue;
            }
            indices.resize(I);
            const UIndirectList<vector> pC
            (
                this->mesh_.Cf().boundaryField()[patchi],
                indices
            );
            UIndirectList<Type> pf(fieldBf[patchi], indices);
            this->getBoundaryField(patchi, indices, pC, pf);
        }
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
Foam::tmp<Foam::Field<Type>> Foam::PointFieldSetType<Type>::getBoundary
(
    const label patchi,
    const GeometricField<Type, pointPatchField, pointMesh>& f
) const
{
    if (isA<valuePointPatchField<Type>>(f.boundaryField()[patchi]))
    {
        return dynamicCast<const valuePointPatchField<Type>>
        (
            f.boundaryField()[patchi]
        );
    }
    return f.boundaryField()[patchi].patchInternalField();
}


template<class Type>
void Foam::PointFieldSetType<Type>::setField()
{
    if (!this->noInternal_)
    {
        const UIndirectList<vector> pts
        (
            this->mesh_.points(),
            this->selectedIndices_
        );
        UIndirectList<Type> fInt(this->fieldPtr_(), this->selectedIndices_);
        this->getInternalField(this->selectedIndices_, pts, fInt);
    }

    labelHashSet points(this->selectedIndices_);

    typename GeometricField<Type, pointPatchField, pointMesh>::
        Boundary& fieldBf = this->fieldPtr_->boundaryFieldRef();

    forAll(this->fieldPtr_->boundaryField(), patchi)
    {
        const polyPatch& p = this->mesh_.boundaryMesh()[patchi];
        if
        (
            this->boundaries_.found(p.name())
         && isA<valuePointPatchField<Type>>(fieldBf[patchi])
        )
        {
            const labelList& meshPoints = p.meshPoints();
            labelList indices(meshPoints.size(), -1);
            labelList gIndices(meshPoints.size(), -1);
            label I = 0;
            forAll(meshPoints, pi)
            {
                if (points.found(meshPoints[pi]))
                {
                    gIndices[I] = meshPoints[pi];
                    indices[I++] = pi;
                }
            }
            indices.resize(I);
            const UIndirectList<vector> pC
            (
                this->mesh_.points(),
                gIndices
            );
            UIndirectList<Type> pf
            (
                dynamicCast<valuePointPatchField<Type>>(fieldBf[patchi]),
                indices
            );
            this->getBoundaryField(patchi, indices, pC, pf);
        }
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

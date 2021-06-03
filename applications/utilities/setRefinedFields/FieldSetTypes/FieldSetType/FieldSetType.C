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

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::FieldSetType<Type, Patch, Mesh>::FieldSetType
(
    const fvMesh& mesh,
    const word& fieldName,
    const labelList& selectedCells,
    Istream& is,
    const bool write
)
:
    mesh_(mesh),
    fieldPtr_(lookupOrRead(fieldName)),
    selectedCells_(selectedCells),
    write_(write),
    good_(fieldPtr_ != nullptr)
{}


template<class Type, template<class> class Patch, class Mesh>
Foam::FieldSetType<Type, Patch, Mesh>::FieldSetType
(
    const fvMesh& mesh,
    const labelList& selectedCells
)
:
    mesh_(mesh),
    fieldPtr_(nullptr),
    selectedCells_(selectedCells),
    write_(false),
    good_(false)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::FieldSetType<Type, Patch, Mesh>::~FieldSetType()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
            new FieldType(fieldHeader, mesh_)
        );
        fPtr->store(fPtr);
        return &mesh_.lookupObjectRef<FieldType>(fieldName);
    }
    return nullptr;
}


template<class Type, template<class> class Patch, class Mesh>
void Foam::FieldSetType<Type, Patch, Mesh>::setField()
{
    typename GeometricField<Type, Patch, Mesh>::
        Boundary& fieldBf = fieldPtr_->boundaryFieldRef();
    forAll(fieldPtr_->boundaryField(), patchi)
    {
        fieldBf[patchi] = fieldBf[patchi].patchInternalField();
    }

    if (write_)
    {
        if (!fieldPtr_->write())
        {
            FatalErrorInFunction
                << "Failed writing field " << fieldPtr_->name() << endl
                << exit(FatalError);
        }
    }
}

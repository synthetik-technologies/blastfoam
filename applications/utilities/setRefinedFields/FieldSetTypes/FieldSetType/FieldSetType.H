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

Class
    Foam::FieldSetType

Description
    Base class used to select how a given field value is computed

SourceFiles
    FieldSetType.C
    FieldSetTypeNew.C

\*---------------------------------------------------------------------------*/

#ifndef FieldSetType_H
#define FieldSetType_H

#include "fvMesh.H"
#include "pointMesh.H"
#include "Time.H"
#include "UautoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class FieldSetType Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class Patch, class Mesh>
class FieldSetType
{
protected:

    typedef Patch<Type> PatchType;
    typedef Mesh MeshType;
    typedef GeometricField<Type, Patch, Mesh> FieldType;

    //- Constant reference to the mesh
    const fvMesh& mesh_;

    //- Reference to the parent dictionary
    const dictionary& dict_;

    //- Name of the field to set
    const word fName_;

    //- Pointer to the field to be set
    UautoPtr<FieldType> fieldPtr_;

    //- Indices to set
    const labelList& selectedIndices_;

    //- Only set patch values
    bool noInternal_;

    //- Evaluate boundaries
    bool evaluateBoundaries_;

    //- Boundaries to set
    wordHashSet boundaries_;

    //- Write after setting fields
    bool write_;

    //- Was the setting successful
    bool good_;

    //- Lookup from the mesh or read the field
    FieldType* lookupOrRead(const word& fieldName) const;

    //- Set the selected cells
    virtual void setField() = 0;

    //- Return the mesh used to construct fields
    const fvMesh& getMesh
    (
        const UautoPtr<GeometricField<Type, fvPatchField, volMesh>>&
    ) const;
    const fvMesh& getMesh
    (
        const UautoPtr<GeometricField<Type, fvsPatchField, surfaceMesh>>&
    ) const;
    const pointMesh& getMesh
    (
        const UautoPtr<GeometricField<Type, pointPatchField, pointMesh>>&
    ) const;

    //- Return the subset of values given indices and positions
    virtual void getInternalField
    (
        const labelList& indices,
        const UIndirectList<vector>& pts,
        UIndirectList<Type>& f
    );

    //- Return the subset of values given indices and positions
    virtual void getBoundaryField
    (
        const label patchi,
        const labelList& indices,
        const UIndirectList<vector>& pts,
        UIndirectList<Type>& f
    );

    //- Return the boundary field value of a field
    virtual tmp<Field<Type>> getBoundary
    (
        const label,
        const FieldType&
    ) const = 0;


public:

    TypeName("FieldSetType");

     // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            FieldSetType,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const word& fieldName,
                const labelList& selectedIndices,
                Istream& is,
                const bool write
            ),
            (mesh, dict, fieldName, selectedIndices, is, write)
        );

    // Constructor
    FieldSetType
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& fieldName,
        const labelList& selectedIndices,
        Istream& is,
        const bool write
    );

    // Constructor
    FieldSetType
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const labelList& selectedIndices
    );


    // Selectors

        static autoPtr<FieldSetType> New
        (
            const word& fieldDesc,
            const word& fieldName,
            const fvMesh& mesh,
            const dictionary& dict,
            const labelList& selectedIndices,
            Istream& is,
            const bool write
        );

    //- Destructor
    virtual ~FieldSetType();

    //- Was the field written successfully
    bool good() const
    {
        return good_;
    }

    //- Return the field to set
    const word& fieldName() const
    {
        return fName_;
    }
};


/*---------------------------------------------------------------------------*\
                          Class VolFieldSetType Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class VolFieldSetType
:
    public FieldSetType<Type, fvPatchField, volMesh>
{
protected:


    //- Set the selected cells
    virtual void setField();

    //- Return the boundary field value of a field
    virtual tmp<Field<Type>> getBoundary
    (
        const label,
        const GeometricField<Type, fvPatchField, volMesh>&
    ) const;

    const fvMesh& getMesh() const
    {
        return this->mesh_;
    }


public:

    // Constructor
    VolFieldSetType
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& fieldName,
        const labelList& selectedIndices,
        Istream& is,
        const bool write
    );

    // Constructor
    VolFieldSetType
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const labelList& selectedIndices
    );

    //- Destructor
    virtual ~VolFieldSetType();
};


/*---------------------------------------------------------------------------*\
                          Class SurfaceFieldSetType Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SurfaceFieldSetType
:
    public FieldSetType<Type, fvsPatchField, surfaceMesh>
{
protected:

    //- Set the selected cells
    virtual void setField();

    //- Return the boundary field value of a field
    virtual tmp<Field<Type>> getBoundary
    (
        const label,
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    ) const;

    const fvMesh& getMesh() const
    {
        return this->mesh_;
    }


public:

    // Constructor
    SurfaceFieldSetType
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& fieldName,
        const labelList& selectedIndices,
        Istream& is,
        const bool write
    );

    // Constructor
    SurfaceFieldSetType
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const labelList& selectedIndices
    );

    //- Destructor
    virtual ~SurfaceFieldSetType();
};


/*---------------------------------------------------------------------------*\
                    Class PointFieldSetType Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PointFieldSetType
:
    public FieldSetType<Type, pointPatchField, pointMesh>
{
protected:

    //- Set the selected cells
    virtual void setField();

    //- Return the boundary field value of a field
    virtual tmp<Field<Type>> getBoundary
    (
        const label,
        const GeometricField<Type, pointPatchField, pointMesh>&
    ) const;

    const pointMesh& getMesh() const
    {
        return pointMesh::New(this->mesh_);
    }

public:

    // Constructor
    PointFieldSetType
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& fieldName,
        const labelList& selectedIndices,
        Istream& is,
        const bool write
    );

    // Constructor
    PointFieldSetType
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const labelList& selectedIndices
    );

    //- Destructor
    virtual ~PointFieldSetType();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFieldSetType(FieldSet)                                             \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(FieldSet, 0);                          \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        FieldSet,                                                              \
        dictionary                                                             \
    );

#define makeFieldSetTypeType(SS, Type, FieldSet)                               \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        FieldSetTypes::SS<Type>,                                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    FieldSet<Type>::adddictionaryConstructorToTable                            \
    <                                                                          \
        FieldSetTypes::SS<Type>                                                \
    >  add##SS##Type##ConstructorToTable_;



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FieldSetType.C"
    #include "FieldSetTypeNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

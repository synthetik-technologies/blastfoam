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
    newFieldSetType.C

\*---------------------------------------------------------------------------*/

#ifndef FieldSetType_H
#define FieldSetType_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvMesh.H"
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

    typedef GeometricField<Type, Patch, Mesh> FieldType;

    //- Constant referenc to the mesh
    const fvMesh& mesh_;

    //- Pointer to the field to be set
    UautoPtr<FieldType> fieldPtr_;

    //- Cells to set
    const labelList& selectedCells_;

    //- Write after setting fields
    bool write_;

    //- Was the setting successful
    bool good_;

    //- Lookup from the mesh or read the field
    FieldType* lookupOrRead(const word& fieldName) const;

    //- Set the selected cells
    virtual void setField();


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
                const word& fieldName,
                const labelList& selectedCells,
                Istream& is,
                const bool write
            ),
            (mesh, fieldName, selectedCells, is, write)
        );

    // Constructor
    FieldSetType
    (
        const fvMesh& mesh,
        const word& fieldName,
        const labelList& selectedCells,
        Istream& is,
        const bool write
    );

    // Constructor
    FieldSetType
    (
        const fvMesh& mesh,
        const labelList& selectedCells
    );


    // Selector
        static autoPtr<FieldSetType> New
        (
            const word& fieldDesc,
            const fvMesh& mesh,
            const labelList& selectedCells,
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
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFieldSetType(Type, FieldSet)                                       \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Type##FieldSet, 0);                    \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        Type##FieldSet,                                                        \
        dictionary                                                             \
    );


#define makeFieldSetTypeType(SS, Type, FieldSet)                               \
                                                                               \
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
    #include "newFieldSetType.C"
#endif

#endif

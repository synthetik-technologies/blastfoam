/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::immersedBoundaryFvPatchField

Description
    This boundary condition is not designed to be evaluated; it is assmued
    that the value is assigned via field assignment, and not via a call to
    e.g. \c updateCoeffs or \c evaluate.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            immersedBoundary;
        value           uniform (0 0 0);    // Required value entry
    }
    \endverbatim

SourceFiles
    immersedBoundaryFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryFvPatchField_H
#define immersedBoundaryFvPatchField_H

#include "volFields.H"
#include "immersedBoundaryObject.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class immersedBoundaryFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class immersedBoundaryFvPatchField
{
protected:

    //- Reference to immersedBoundaryPolyPatch
    const immersedBoundaryObject& ibm_;

    //- Reference to field
    GeometricField<Type, fvPatchField, volMesh>& field_;

    //- Values on patch
    mutable Field<Type> values_;

    //- Set the internal cells to a given value
    void setInternal(const Type& val);

    //- Set the average value for all internal cells
    void setInternalAvg();

    //- Set the average value for based on the face values
    void setPatchAvg();

    //- Smooth the internal cells
    void setSmoothInternal();


public:

    //- Runtime type information
    TypeName("immersedBoundaryFvPatchField");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            immersedBoundaryFvPatchField,
            dictionary,
            (
                GeometricField<Type, fvPatchField, volMesh>& f,
                const dictionary& dict,
                const immersedBoundaryObject& ibo
            ),
            (f, dict, ibo)
        );

    virtual ~immersedBoundaryFvPatchField()
    {}

    // Constructors

        //- Construct by mapping given patchField<Type> onto a new patch
        immersedBoundaryFvPatchField
        (
            GeometricField<Type, fvPatchField, volMesh>&,
            const dictionary&,
            const immersedBoundaryObject&
        );


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static autoPtr<immersedBoundaryFvPatchField<Type>> New
        (
            GeometricField<Type, fvPatchField, volMesh>&,
            const dictionary&,
            const immersedBoundaryObject&
        );


    // Member Functions


        // Evaluation functions

            //- Update size and values
            virtual void updateCoeffs() const = 0;

            //- Add forcing to the given field
            virtual void addForcing
            (
                Field<Type>& f,
                const Field<scalar>& alphaRho,
                const Field<Type>& old,
                const Field<Type>& RHS,
                const scalar& dt
            ) const;


            //- Add forcing to the given field
            virtual void setValues() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "immersedBoundaryFvPatchField.C"
#endif

#define makeImmersedFvPatchField(PatchTypeField)                              \
    defineNamedTemplateTypeNameAndDebug(PatchTypeField, 0);                   \
    defineTemplateRunTimeSelectionTable(PatchTypeField, dictionary);


#define addToImmersedPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField) \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary                                                             \
    );

// For non-templated patch fields
#define makeImmersedPatchTypeField(PatchTypeField, typePatchTypeField)         \
    defineTypeNameAndDebug(typePatchTypeField, 0);                             \
    addToImmersedPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)

// For templated patch fields
#define makeTemplateImmersedPatchTypeField(fieldType, type)                    \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, ImmersedBoundaryFvPatch, CAPITALIZE(fieldType), Field),     \
        0                                                                      \
    );                                                                         \
    addToImmersedPatchFieldRunTimeSelection                                    \
    (                                                                          \
        CAT3(immersedBoundary, CAPITALIZE(fieldType), PatchField),             \
        CAT4(type, ImmersedBoundaryFvPatch, CAPITALIZE(fieldType), Field)      \
    )

#define makeImmersedPatchFields(type)                                          \
    FOR_ALL_FIELD_TYPES(makeTemplateImmersedPatchTypeField, type)

#define makeImmersedPatchFieldTypeName(fieldType, type)                        \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, ImmersedBoundaryFvPatch, CAPITALIZE(fieldType), Field),     \
        0                                                                      \
    );

#define makeImmersedPatchFieldTypeNames(type)                                  \
    FOR_ALL_FIELD_TYPES(makeImmersedPatchFieldTypeName, type)

#define makeImmersedPatchTypeFieldTypedef(fieldType, type)                     \
    typedef type##ImmersedBoundaryFvPatchField<fieldType>                      \
        CAT4(type, ImmersedBoundaryFvPatch, CAPITALIZE(fieldType), Field);

#define makeImmersedPatchTypeFieldTypedefs(type)                               \
    FOR_ALL_FIELD_TYPES(makeImmersedPatchTypeFieldTypedef, type)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

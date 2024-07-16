/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Class
    Foam::fePatchField

Description
    Abstract base class for fe-mesh patch fields.

    The base-field does not store values as they are part of the
    "internal field".  There are derived classes to store constraint values
    e.g. fixedValueFePatchField derived from the generic
    valueFePatchField which ensures the values in the "internal field"
    are reset to the fixed-values by applying the stored values.

SourceFiles
    fePatchField.C
    newFePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fePatchField_H
#define fePatchField_H

#include "fePatch.H"
#include "DimensionedField.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class objectRegistry;
class dictionary;
class fePatchFieldMapper;
class feMesh;

// Forward declaration of friend functions and operators

template<class Type>
class fePatchField;

template<class Type>
class calculatedFePatchField;

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const fePatchField<Type>&
);


/*---------------------------------------------------------------------------*\
                       Class fePatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fePatchField
{
    // Private Data

        //- Reference to patch
        const fePatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, feMesh>& internalField_;

        //- Update index used so that updateCoeffs is called only once during
        //  the construction of the matrix
        bool updated_;


public:

    typedef Type value_type;
    typedef fePatch Patch;
    typedef calculatedFePatchField<Type> Calculated;


    //- Runtime type information
    TypeName("fePatchField");

    //- Debug switch to disallow the use of genericFePatchField
    static int disallowGenericFePatchField;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            fePatchField,
            fePatch,
            (
                const fePatch& p,
                const DimensionedField<Type, feMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            fePatchField,
            patchMapper,
            (
                const fePatchField<Type>& ptf,
                const fePatch& p,
                const DimensionedField<Type, feMesh>& iF,
                const fePatchFieldMapper& m
            ),
            (dynamic_cast<const fePatchFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            fePatchField,
            dictionary,
            (
                const fePatch& p,
                const DimensionedField<Type, feMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        fePatchField
        (
            const fePatch&,
            const DimensionedField<Type, feMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fePatchField
        (
            const fePatch&,
            const DimensionedField<Type, feMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        fePatchField
        (
            const fePatchField<Type>&,
            const fePatch&,
            const DimensionedField<Type, feMesh>&,
            const fePatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        fePatchField(const fePatchField<Type>&) = delete;

        //- Disallow clone without setting internal field reference
        autoPtr<fePatchField<Type>> clone() const
        {
            NotImplemented;
            return tmp<fePatchField<Type>>(nullptr);
        }

        //- Construct as copy setting internal field reference
        fePatchField
        (
            const fePatchField<Type>&,
            const DimensionedField<Type, feMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fePatchField<Type>> clone
        (
            const DimensionedField<Type, feMesh>& iF
        ) const = 0;


    // Selectors

        //- Return a feer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static autoPtr<fePatchField<Type>> New
        (
            const word&,
            const fePatch&,
            const DimensionedField<Type, feMesh>&
        );

        //- Return a feer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values).
        //  Allows override of constraint type
        static autoPtr<fePatchField<Type>> New
        (
            const word&,
            const word& actualPatchType,
            const fePatch&,
            const DimensionedField<Type, feMesh>&
        );

        //- Return a feer to a new patchField created on freestore from
        //  a given fePatchField mapped onto a new patch
        static autoPtr<fePatchField<Type>> New
        (
            const fePatchField<Type>&,
            const fePatch&,
            const DimensionedField<Type, feMesh>&,
            const fePatchFieldMapper&
        );

        //- Return a feer to a new patchField created on freestore
        //  from dictionary
        static autoPtr<fePatchField<Type>> New
        (
            const fePatch&,
            const DimensionedField<Type, feMesh>&,
            const dictionary&
        );

        //- Return a feer to a new calculatedFePatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static autoPtr<fePatchField<Type>>
        NewCalculatedType
        (
            const fePatchField<Type2>&
        );


    //- Destructor
    virtual ~fePatchField<Type>()
    {}


    // Member Functions

        // Attributes

            //- Return the type of the calculated form of fePatchField
            static const word& calculatedType();

            //- Return size
            label size() const
            {
                return patch().size();
            }

            //- Return true if this patch field fixes a value
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return false;
            }

            //- Return true if this overrides the underlying constraint type
            bool overridesConstraint() const
            {
                if (type() == patch_.type())
                {
                    return false;
                }

                typename fePatchConstructorTable::iterator patchTypeCstrIter
                    = fePatchConstructorTablePtr_->find(patch_.type());

                return
                    patchTypeCstrIter
                 != fePatchConstructorTablePtr_->end();
            }


        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Return patch
            const fePatch& patch() const
            {
                return patch_;
            }

            //- Return dimensioned internal field reference
            const DimensionedField<Type, feMesh>& internalField() const
            {
                return internalField_;
            }

            //- Return internal field reference
            const Field<Type>& primitiveField() const
            {
                return internalField_;
            }

            //- Return true if the boundary condition has already been updated
            bool updated() const
            {
                return updated_;
            }

            //- Return field created from appropriate internal field values
            tmp<Field<Type>> patchInternalField() const;

            //- Return field created from appropriate internal field values
            //  given internal field reference
            template<class Type1>
            tmp<Field<Type1>> patchInternalField
            (
                const Field<Type1>& iF
            ) const;

            //- Return field created from selected internal field values
            //  given internal field reference
            template<class Type1>
            tmp<Field<Type1>> patchInternalField
            (
                const Field<Type1>& iF,
                const labelList& meshFes
            ) const;

            //- Given the internal field and a patch field,
            //  add the patch field to the internal field
            template<class Type1>
            void addToInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF
            ) const;

            //- Given the internal field and a patch field,
            //  add selected elements of the patch field to the internal field
            template<class Type1>
            void addToInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF,
                const labelList& fes
            ) const;

            //- Given the internal field and a patch field,
            //  set the patch field in the internal field
            template<class Type1>
            void setInInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF,
                const labelList& meshFes
            ) const;

            //- Given the internal field and a patch field,
            //  set the patch field in the internal field
            template<class Type1>
            void setInInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF
            ) const;


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fePatchFieldMapper&)
            {}

            //- Reverse map the given fePatchField onto this fePatchField
            //  Used to reconstruct fields
            virtual void rmap(const fePatchField<Type>&, const labelList&)
            {}


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs()
            {
                updated_ = true;
            }

            //- Initialise evaluation of the patch field (do nothing)
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            )
            {}

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;


    // Member Operators

        virtual void operator=(const fePatchField<Type>&){}
        virtual void operator=(fePatchField<Type>&&){}
        virtual void operator+=(const fePatchField<Type>&){}
        virtual void operator-=(const fePatchField<Type>&){}
        virtual void operator*=(const fePatchField<scalar>&){}
        virtual void operator/=(const fePatchField<scalar>&){}

        virtual void operator=(const Field<Type>&){}
        virtual void operator+=(const Field<Type>&){}
        virtual void operator-=(const Field<Type>&){}

        virtual void operator*=(const Field<scalar>&){}
        virtual void operator/=(const Field<scalar>&){}

        virtual void operator=(const Type&){}
        virtual void operator+=(const Type&){}
        virtual void operator-=(const Type&){}
        virtual void operator*=(const scalar){}
        virtual void operator/=(const scalar){}


        // Force an assignment irrespective of form of patch
        // By generic these do nothing unless the patch actually has boundary
        // values

            virtual void operator==(const fePatchField<Type>&){}
            virtual void operator==(const Field<Type>&){}
            virtual void operator==(const Type&){}


    // Ostream operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const fePatchField<Type>&
        );
};


// This function is added as a hack to enable simple backward compatibility
// with versions using referenceLevel in GeometicField
template<class Type>
const fePatchField<Type>& operator+
(
    const fePatchField<Type>& ppf,
    const Type&
)
{
    return ppf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fePatchFieldFunctions.H"

#ifdef NoRepository
    #include "fePatchField.C"
    #include "calculatedFePatchField.H"
#endif


#define addToFePatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField) \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        fePatch                                                             \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patchMapper                                                            \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary                                                             \
    );


// for non-templated patch fields
#define makeFePatchTypeField(PatchTypeField,typePatchTypeField)             \
    defineTypeNameAndDebug(typePatchTypeField, 0);                             \
    addToFePatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)


// for templated patch fields
#define makeTemplateFePatchTypeField(PatchTypeField, typePatchTypeField)    \
    defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);                \
    addToFePatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)


#define makeFePatchFields(type)                                             \
    makeTemplateFePatchTypeField                                            \
    (                                                                          \
        fePatchScalarField,                                                 \
        type##FePatchScalarField                                            \
    );                                                                         \
    makeTemplateFePatchTypeField                                            \
    (                                                                          \
        fePatchVectorField,                                                 \
        type##FePatchVectorField                                            \
    );                                                                         \
    makeTemplateFePatchTypeField                                            \
    (                                                                          \
        fePatchSphericalTensorField,                                        \
        type##FePatchSphericalTensorField                                   \
    );                                                                         \
    makeTemplateFePatchTypeField                                            \
    (                                                                          \
        fePatchSymmTensorField,                                             \
        type##FePatchSymmTensorField                                        \
    );                                                                         \
    makeTemplateFePatchTypeField                                            \
    (                                                                          \
        fePatchTensorField,                                                 \
        type##FePatchTensorField                                            \
);


#define makeFePatchFieldsTypeName(type)                                     \
    defineNamedTemplateTypeNameAndDebug(type##FePatchScalarField, 0);       \
    defineNamedTemplateTypeNameAndDebug(type##FePatchVectorField, 0);       \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        type##FePatchSphericalTensorField, 0                                \
    );                                                                         \
    defineNamedTemplateTypeNameAndDebug(type##FePatchSymmTensorField, 0);   \
    defineNamedTemplateTypeNameAndDebug(type##FePatchTensorField, 0)


#define makeFePatchFieldTypedefs(type)                                      \
    typedef type##FePatchField<scalar> type##FePatchScalarField;         \
    typedef type##FePatchField<vector> type##FePatchVectorField;         \
    typedef type##FePatchField<sphericalTensor>                             \
        type##FePatchSphericalTensorField;                                  \
    typedef type##FePatchField<symmTensor> type##FePatchSymmTensorField;\
    typedef type##FePatchField<tensor> type##FePatchTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

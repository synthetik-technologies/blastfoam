/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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
    Foam::PatchFieldSelector

Description
    Class to select the field from a given patchField

SourceFiles
    PatchFieldSelector.C
    PatchFieldSelector.C

\*---------------------------------------------------------------------------*/

#ifndef PatchFieldSelector_H
#define PatchFieldSelector_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PtrList.H"
#include "fvPatchField.H"
#include "fvsPatchField.H"
#include "pointPatchField.H"

#include "calculatedFvPatchField.H"
#include "fixedValueFvPatchField.H"
#include "fixedGradientFvPatchField.H"
#include "mixedFvPatchField.H"

#include "calculatedFvsPatchField.H"
#include "fixedValueFvsPatchField.H"

#include "valuePointPatchField.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class patchFieldSelectorBase Declaration
\*---------------------------------------------------------------------------*/

class patchFieldSelectorBase
{
public:
    ClassName("patchFieldSelector");

    patchFieldSelectorBase()
    {}
};

/*---------------------------------------------------------------------------*\
                    Class PatchFieldSelector Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PatchFieldSelector
:
    public patchFieldSelectorBase
{
public:

    // Constructor
    PatchFieldSelector()
    {}

    // Selectors

        //- Selector for fvPatchField
        static PtrList<PatchFieldSelector<Type>> New
        (
            fvPatchField<Type>& pfield,
            const dictionary& dict
        );

        //- Selector for fvsPatchField
        static PtrList<PatchFieldSelector<Type>> New
        (
            fvsPatchField<Type>& pfield,
            const dictionary& dict
        );

        //- Selector for pointPatchField
        static PtrList<PatchFieldSelector<Type>> New
        (
            pointPatchField<Type>& pfield,
            const dictionary& dict
        );


    //- Destructor
    virtual ~PatchFieldSelector()
    {}


    // Member functions

        //- Return the patch index
        virtual label index() const = 0;

        //- Return the type
        virtual word type() const
        {
            return word::null;
        }

        //- Return the field to relax given a fvPatchField
        virtual const Field<Type>& relaxField
        (
            const fvPatchField<Type>& f
        ) const
        {
            NotImplemented;
            return f;
        }

        //- Return the field to relax given a fvsPatchField
        virtual const Field<Type>& relaxField
        (
            const fvsPatchField<Type>& f
        ) const
        {
            NotImplemented;
            return f;
        }

        //- Return the field to relax given a pointPatchField
        virtual const Field<Type>& relaxField
        (
            const pointPatchField<Type>& f
        ) const
        {
            NotImplemented;
            return dynamicCast<const Field<Type>>(f);
        }

        //- Return the field to relax
        virtual const Field<Type>& field() const = 0;

        //- Access the field to relax
        virtual Field<Type>& field() = 0;

        //- Store the previous iteration for the give field
        virtual void storePrevIter()
        {}

        //- Correct the underlying field
        virtual void correct()
        {}
};


/*---------------------------------------------------------------------------*\
                Class ValueFvPatchFieldSelector Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ValueFvPatchFieldSelector
:
    public PatchFieldSelector<Type>
{
// Private member data

    //- Reference to the fvPatchField
    fvPatchField<Type>& pfield_;


public:

    // Constructor
    ValueFvPatchFieldSelector(fvPatchField<Type>& pfield)
    :
        pfield_(pfield)
    {}


    //- Destructor
    virtual ~ValueFvPatchFieldSelector()
    {}


    // Member functions

        //- Return the patch index
        inline label index() const
        {
            return pfield_.patch().index();
        }

        //- Return the field to relax given an fvPatchField
        using PatchFieldSelector<Type>::relaxField;
        inline const Field<Type>& relaxField
        (
            const fvPatchField<Type>& f
        ) const
        {
            return f;
        }

        //- Return the field to relax
        inline const Field<Type>& field() const
        {
            return pfield_;
        }

        //- Access the field to relax
        inline Field<Type>& field()
        {
            return pfield_;
        }
};


/*---------------------------------------------------------------------------*\
            Class FixedGradientFvPatchFieldSelector Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class FixedGradientFvPatchFieldSelector
:
    public PatchFieldSelector<Type>
{
// Private member data

    //- Reference to the fvPatchField
    fixedGradientFvPatchField<Type>& pfield_;


public:

    // Constructor
    FixedGradientFvPatchFieldSelector(fvPatchField<Type>& pfield)
    :
        pfield_(dynamicCast<fixedGradientFvPatchField<Type>>(pfield))
    {}


    //- Destructor
    virtual ~FixedGradientFvPatchFieldSelector()
    {}


    // Member functions

        //- Return the patch index
        inline label index() const
        {
            return pfield_.patch().index();
        }

        //- Return the type
        virtual word type() const
        {
            return "(grad)";
        }

        //- Return the field to relax given another patch field
        using PatchFieldSelector<Type>::relaxField;
        inline const Field<Type>& relaxField
        (
            const fvPatchField<Type>& f
        ) const
        {
            return
                dynamicCast<const fixedGradientFvPatchField<Type>>
                (
                    f
                ).gradient();
        }

        //- Return the field to relax
        inline const Field<Type>& field() const
        {
            return pfield_.gradient();
        }

        //- Access the field to relax
        inline Field<Type>& field()
        {
            return pfield_.gradient();
        }

        //- Store the previous iteration for the give field
        virtual void storePrevIter()
        {
            const fvPatchField<Type>& pfieldPrev =
                dynamicCast<const GeometricField<Type, fvPatchField, volMesh>>
                (
                    pfield_.internalField()
                ).prevIter().boundaryField()[index()];
            dynamicCast<fixedGradientFvPatchField<Type>>
            (
                const_cast<fvPatchField<Type>&>(pfieldPrev)
            ).gradient() = field();
        }

        //- Correct the underlying field
        virtual void correct()
        {
            pfield_.fvPatchField<Type>::updateCoeffs();
            pfield_.evaluate();
        }
};


/*---------------------------------------------------------------------------*\
                Class MixedFvPatchFieldSelector Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MixedFvPatchFieldSelector
:
    public PatchFieldSelector<Type>
{
public:
// Public enumeration
    enum FieldType
    {
        REF_VALUE,
        REF_GRADIENT,
        VALUE_FRACTION
    };


private:
// Private member data

    //- Reference to the fvPatchField
    mixedFvPatchField<Type>& pfield_;

    //- Field to return
    const FieldType type_;


public:

    // Constructor
    MixedFvPatchFieldSelector
    (
        fvPatchField<Type>& pfield,
        const FieldType type
    )
    :
        pfield_(dynamicCast<mixedFvPatchField<Type>>(pfield)),
        type_(type)
    {}


    //- Destructor
    virtual ~MixedFvPatchFieldSelector()
    {}


    // Member functions

        //- Return the patch index
        inline label index() const
        {
            return pfield_.patch().index();
        }

        //- Return the type
        virtual word type() const
        {
            switch (type_)
            {
                case REF_VALUE:
                    return "(refValue)";
                case REF_GRADIENT:
                    return "(refGrad)";
                case VALUE_FRACTION:
                    NotImplemented;
                    return "(valueFraction)";//mpf.valueFraction();
            }
            return word::null;
        }

        //- Return the field to relax given an fvPatchField
        using PatchFieldSelector<Type>::relaxField;
        inline const Field<Type>& relaxField
        (
            const fvPatchField<Type>& f
        ) const
        {
            const mixedFvPatchField<Type>& mpf =
                dynamicCast<const mixedFvPatchField<Type>>(f);
            switch (type_)
            {
                case REF_VALUE:
                    return mpf.refValue();
                case REF_GRADIENT:
                    return mpf.refGrad();
                case VALUE_FRACTION:
                    NotImplemented;
                    return f;//mpf.valueFraction();
            }
            return f;
        }

        //- Return the field to relax
        inline const Field<Type>& field() const
        {
            switch (type_)
            {
                case REF_VALUE:
                    return pfield_.refValue();
                case REF_GRADIENT:
                    return pfield_.refGrad();
                case VALUE_FRACTION:
                    NotImplemented;
                    return pfield_;//pfield_.valueFraction();
            }
            return pfield_;
        }

        //- Access the field to relax
        inline Field<Type>& field()
        {
            switch (type_)
            {
                case REF_VALUE:
                    return pfield_.refValue();
                case REF_GRADIENT:
                    return pfield_.refGrad();
                case VALUE_FRACTION:
                    NotImplemented;
                    return pfield_;//pfield_.valueFraction();
            }
            return pfield_;
        }

        //- Store the previous iteration for the give field
        virtual void storePrevIter()
        {
            const fvPatchField<Type>& pfieldPrev =
                dynamicCast<const GeometricField<Type, fvPatchField, volMesh>>
                (
                    pfield_.internalField()
                ).prevIter().boundaryField()[index()];
            mixedFvPatchField<Type>& pmixedPrev =
                dynamicCast<mixedFvPatchField<Type>>
                (
                    const_cast<fvPatchField<Type>&>(pfieldPrev)
                );
            switch (type_)
            {
                case REF_VALUE:
                    pmixedPrev.refValue() = pfield_.refValue();
                    return;
                case REF_GRADIENT:
                    pmixedPrev.refGrad() = pfield_.refGrad();
                    return;
                case VALUE_FRACTION:
                    NotImplemented;
                    pmixedPrev.valueFraction() = pfield_.valueFraction();
                    return;
            }
        }

        //- Correct the underlying field
        virtual void correct()
        {
            pfield_.fvPatchField<Type>::updateCoeffs();
            pfield_.evaluate();
        }
};


/*---------------------------------------------------------------------------*\
                Class ValueFvPatchFieldSelector Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ValueFvsPatchFieldSelector
:
    public PatchFieldSelector<Type>
{
// Private member data

    //- Reference to the fvPatchField
    fvsPatchField<Type>& pfield_;


public:

    // Constructor
    ValueFvsPatchFieldSelector(fvsPatchField<Type>& pfield)
    :
        pfield_(pfield)
    {}


    //- Destructor
    virtual ~ValueFvsPatchFieldSelector()
    {}


    // Member functions

        //- Return the patch index
        inline label index() const
        {
            return pfield_.patch().index();
        }

        //- Return the field to relax given an fvPatchField
        using PatchFieldSelector<Type>::relaxField;
        inline const Field<Type>& relaxField
        (
            const fvsPatchField<Type>& f
        ) const
        {
            return f;
        }

        //- Return the field to relax
        inline const Field<Type>& field() const
        {
            return pfield_;
        }

        //- Access the field to relax
        inline Field<Type>& field()
        {
            return pfield_;
        }
};


/*---------------------------------------------------------------------------*\
                Class ValuePointPatchFieldSelector Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ValuePointPatchFieldSelector
:
    public PatchFieldSelector<Type>
{
// Private member data

    //- Reference to the fvPatchField
    valuePointPatchField<Type>& pfield_;


public:

    // Constructor
    ValuePointPatchFieldSelector(pointPatchField<Type>& pfield)
    :
        pfield_(dynamicCast<valuePointPatchField<Type>>(pfield))
    {}


    //- Destructor
    virtual ~ValuePointPatchFieldSelector()
    {}


    // Member functions

        //- Return the patch index
        inline label index() const
        {
            return pfield_.patch().index();
        }

        //- Return the field to relax given an fvPatchField
        using PatchFieldSelector<Type>::relaxField;
        inline const Field<Type>& relaxField
        (
            const pointPatchField<Type>& f
        ) const
        {
            return dynamicCast<const valuePointPatchField<Type>>(f);
        }

        //- Return the field to relax
        inline const Field<Type>& field() const
        {
            return pfield_;
        }

        //- Access the field to relax
        inline Field<Type>& field()
        {
            return pfield_;
        }

        //- Correct the underlying field
        virtual void correct()
        {
            pfield_.setInInternalField
            (
                const_cast<Field<Type>&>(pfield_.primitiveField()),
                pfield_
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchFieldSelector.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

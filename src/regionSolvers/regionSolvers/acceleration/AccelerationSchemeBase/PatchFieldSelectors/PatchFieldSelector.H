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
    //- Index of relaxation patch
    const label index_;

public:

    // Constructor
    PatchFieldSelector(const label index)
    :
        index_(index)
    {}

    // Selectors

        //- Selector for fvPatchField
        static PtrList<PatchFieldSelector<Type>> New
        (
            const fvPatchField<Type>& pfield,
            const dictionary& dict
        );

        //- Selector for fvsPatchField
        static PtrList<PatchFieldSelector<Type>> New
        (
            const fvsPatchField<Type>& pfield,
            const dictionary& dict
        );

        //- Selector for pointPatchField
        static PtrList<PatchFieldSelector<Type>> New
        (
            const pointPatchField<Type>& pfield,
            const dictionary& dict
        );


    //- Destructor
    virtual ~PatchFieldSelector()
    {}


    // Member functions

        //- Return the patch index
        inline label index() const
        {
            return index_;
        }

        //- Return the type
        virtual word type() const
        {
            return word::null;
        }

        // Return relaxation field

            //- Return the field to relax given a fvPatchField
            virtual const Field<Type>& relaxField
            (
                const fvPatchField<Type>& pfield
            ) const
            {
                NotImplemented;
                return pfield;
            }

            //- Access the field to relax given a fvPatchField
            virtual Field<Type>& relaxField
            (
                fvPatchField<Type>& pfield
            ) const
            {
                NotImplemented;
                return pfield;
            }

            //- Return the field to relax given a fvsPatchField
            virtual const Field<Type>& relaxField
            (
                const fvsPatchField<Type>& pfield
            ) const
            {
                NotImplemented;
                return pfield;
            }

            //- Access the field to relax given a fvsPatchField
            virtual Field<Type>& relaxField
            (
                fvsPatchField<Type>& pfield
            ) const
            {
                NotImplemented;
                return pfield;
            }

            //- Return the field to relax given a pointPatchField
            virtual const Field<Type>& relaxField
            (
                const pointPatchField<Type>& pfield
            ) const
            {
                NotImplemented;
                return dynamicCast<const Field<Type>>(pfield);
            }

            //- Access the field to relax given a pointPatchField
            virtual Field<Type>& relaxField
            (
                pointPatchField<Type>& pfield
            ) const
            {
                NotImplemented;
                return dynamicCast<Field<Type>>(pfield);
            }

        // Store the previous iteration

            //- Store the previous iteration for the given fvPatchField
            virtual void storePrevIter
            (
                const fvPatchField<Type>& pfield,
                fvPatchField<Type>& pfieldPrev
            )
            {
                this->relaxField(pfieldPrev) = this->relaxField(pfield);
            }

            //- Store the previous iteration for the given fvPatchField
            virtual void storePrevIter
            (
                const fvsPatchField<Type>& pfield,
                fvsPatchField<Type>& pfieldPrev
            )
            {
                this->relaxField(pfieldPrev) = this->relaxField(pfield);
            }

            //- Store the previous iteration for the given fvPatchField
            virtual void storePrevIter
            (
                const pointPatchField<Type>& pfield,
                pointPatchField<Type>& pfieldPrev
            )
            {
                this->relaxField(pfieldPrev) = this->relaxField(pfield);
            }


        // Correct patch fields (evaluate)

            //- Correct the fvPatchField
            virtual void correct(fvPatchField<Type>& pfield)
            {
                pfield.fvPatchField<Type>::updateCoeffs();
                pfield.evaluate();
            }

            //- Correct the fvsPatchField
            virtual void correct(fvsPatchField<Type>& pfield)
            {}

            //- Correct the pointPatchField
            virtual void correct(pointPatchField<Type>& pfield)
            {
                valuePointPatchField<Type>& vpfield =
                    dynamicCast<valuePointPatchField<Type>>(pfield);
                vpfield.setInternalField
                (
                    const_cast<Field<Type>&>(vpfield.primitiveField()),
                    vpfield
                );
            }
};


/*---------------------------------------------------------------------------*\
                Class ValueFvPatchFieldSelector Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ValueFvPatchFieldSelector
:
    public PatchFieldSelector<Type>
{
public:

    // Constructor
    ValueFvPatchFieldSelector(const fvPatchField<Type>& pfield)
    :
        PatchFieldSelector<Type>(pfield.patch().index())
    {}


    //- Destructor
    virtual ~ValueFvPatchFieldSelector()
    {}


    // Member functions

        //- Return the field to relax given an fvPatchField
        inline const Field<Type>& relaxField
        (
            const fvPatchField<Type>& f
        ) const
        {
            return f;
        }

        //- Access the field to relax given an fvPatchField
        inline Field<Type>& relaxField
        (
            fvPatchField<Type>& f
        ) const
        {
            return f;
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

public:

    // Constructor
    FixedGradientFvPatchFieldSelector(const fvPatchField<Type>& pfield)
    :
        PatchFieldSelector<Type>(pfield.patch().index())
    {}


    //- Destructor
    virtual ~FixedGradientFvPatchFieldSelector()
    {}


    // Member functions

        //- Return the type
        virtual word type() const
        {
            return "(grad)";
        }

        using PatchFieldSelector<Type>::relaxField;

        //- Return the field to relax given another patch field
        inline const Field<Type>& relaxField
        (
            const fvPatchField<Type>& pfield
        ) const
        {
            return
                dynamicCast<const fixedGradientFvPatchField<Type>>
                (
                    pfield
                ).gradient();
        }

        //- Access the field to relax given another patch field
        inline Field<Type>& relaxField(fvPatchField<Type>& pfield) const
        {
            return
                dynamicCast<fixedGradientFvPatchField<Type>>
                (
                    pfield
                ).gradient();
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

    //- Field to return
    const FieldType type_;


public:

    // Constructor
    MixedFvPatchFieldSelector
    (
        const fvPatchField<Type>& pfield,
        const FieldType type
    )
    :
        PatchFieldSelector<Type>(pfield.patch().index()),
        type_(type)
    {}


    //- Destructor
    virtual ~MixedFvPatchFieldSelector()
    {}


    // Member functions

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

        using PatchFieldSelector<Type>::relaxField;

        //- Return the field to relax given an fvPatchField
        inline const Field<Type>& relaxField
        (
            const fvPatchField<Type>& pfield
        ) const
        {
            const mixedFvPatchField<Type>& mixedpfield =
                dynamicCast<const mixedFvPatchField<Type>>(pfield);
            switch (type_)
            {
                case REF_VALUE:
                    return mixedpfield.refValue();
                case REF_GRADIENT:
                    return mixedpfield.refGrad();
                case VALUE_FRACTION:
                    NotImplemented;
                    return pfield;//mixedpfield.valueFraction();
            }
            return pfield;
        }

        //- Access the field to relax given an fvPatchField
        inline Field<Type>& relaxField(fvPatchField<Type>& pfield) const
        {
            mixedFvPatchField<Type>& mixedpfield =
                dynamicCast<mixedFvPatchField<Type>>(pfield);
            switch (type_)
            {
                case REF_VALUE:
                    return mixedpfield.refValue();
                case REF_GRADIENT:
                    return mixedpfield.refGrad();
                case VALUE_FRACTION:
                    NotImplemented;
                    return pfield;//mixedpfield.valueFraction();
            }
            return pfield;
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
public:

    // Constructor
    ValueFvsPatchFieldSelector(const fvsPatchField<Type>& pfield)
    :
        PatchFieldSelector<Type>(pfield.patch().index())
    {}


    //- Destructor
    virtual ~ValueFvsPatchFieldSelector()
    {}


    // Member functions

        using PatchFieldSelector<Type>::relaxField;

        //- Return the field to relax given an fvPatchField
        inline const Field<Type>& relaxField
        (
            const fvsPatchField<Type>& f
        ) const
        {
            return f;
        }

        //- Access the field to relax given an fvPatchField
        inline Field<Type>& relaxField(fvsPatchField<Type>& f) const
        {
            return f;
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
public:

    // Constructor
    ValuePointPatchFieldSelector(const pointPatchField<Type>& pfield)
    :
        PatchFieldSelector<Type>(pfield.patch().index())
    {}


    //- Destructor
    virtual ~ValuePointPatchFieldSelector()
    {}


    // Member functions

        using PatchFieldSelector<Type>::relaxField;

        //- Return the field to relax given an fvPatchField
        inline const Field<Type>& relaxField
        (
            const pointPatchField<Type>& f
        ) const
        {
            return dynamicCast<const valuePointPatchField<Type>>(f);
        }

        //- Access the field to relax given an fvPatchField
        inline Field<Type>& relaxField(pointPatchField<Type>& f) const
        {
            return dynamicCast<valuePointPatchField<Type>>(f);
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

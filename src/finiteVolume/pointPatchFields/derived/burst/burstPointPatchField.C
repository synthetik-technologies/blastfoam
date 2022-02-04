/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "burstPointPatchField.H"
#include "calculatedPointPatchField.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::burstPointPatchField<Type>::burstPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    valuePointPatchField<Type>(p, iF),
    burstPointPatchFieldBase
    (
        dynamicCast<const facePointPatch>(this->patch()).patch()
    ),
    burstPointPatchField_
    (
        new calculatedPointPatchField<Type>
        (
            p,
            iF
        )
    ),
    intactPointPatchField_
    (
        new calculatedPointPatchField<Type>
        (
            p,
            iF
        )
    )
{}


template<class Type>
Foam::burstPointPatchField<Type>::burstPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    valuePointPatchField<Type>(p, iF, dict),
    burstPointPatchFieldBase
    (
        dynamicCast<const facePointPatch>(this->patch()).patch()
    ),
    burstPointPatchField_(),
    intactPointPatchField_()
{
    if (!isType<burstPointPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "patch " << this->patch().index() << " not burst type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }

    // Create a new patch dictionary and replace the type with the intactType
    {
        dictionary burstDict(dict.parent(), dict.subDict("burstPatch"));
        if (!burstDict.found("patchType"))
        {
            burstDict.add("patchType", typeName);
        }
        burstPointPatchField_.set
        (
            pointPatchField<Type>::New
            (
                p,
                iF,
                burstDict
            ).ptr()
        );
    }
    {
        dictionary intactDict(dict.parent(), dict.subDict("intactPatch"));
        if (!intactDict.found("patchType"))
        {
            intactDict.add("patchType", typeName);
        }
        intactPointPatchField_.set
        (
            pointPatchField<Type>::New
            (
                p,
                iF,
                intactDict
            ).ptr()
        );
    }

}


template<class Type>
Foam::burstPointPatchField<Type>::burstPointPatchField
(
    const burstPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    valuePointPatchField<Type>(ptf, p, iF, mapper),
    burstPointPatchFieldBase
    (
        dynamicCast<const facePointPatch>(this->patch()).patch()
    ),
    intactPointPatchField_(ptf.intactPointPatchField_->clone(iF).ptr())
{
    if (!isType<burstPointPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::burstPointPatchField<Type>::burstPointPatchField
(
    const burstPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    valuePointPatchField<Type>(ptf, iF),
    burstPointPatchFieldBase
    (
        dynamicCast<const facePointPatch>(this->patch()).patch()
    ),
    intactPointPatchField_(ptf.intactPointPatchField_->clone(iF).ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::burstPointPatchField<Type>::autoMap
(
    const pointPatchFieldMapper& m
)
{
    pointPatchField<Type>::autoMap(m);
    intactPointPatchField_->autoMap(m);
    burstPointPatchField_->autoMap(m);
}


template<class Type>
void Foam::burstPointPatchField<Type>::rmap
(
    const pointPatchField<Type>& ptf,
    const labelList& addr
)
{
    pointPatchField<Type>::rmap(ptf, addr);

    const burstPointPatchField<Type>& bpf =
        refCast<const burstPointPatchField<Type>>(ptf);
    intactPointPatchField_->rmap(bpf.intactPointPatchField_(), addr);
    burstPointPatchField_->rmap(bpf.burstPointPatchField_(), addr);
}


template<class Type>
void Foam::burstPointPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    intactPointPatchField_->updateCoeffs();
    burstPointPatchField_->updateCoeffs();
}


template<class Type>
void Foam::burstPointPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Save the origin internal values;
    Field<Type> pIf(this->patchInternalField());

    // Get internal field to insert values into
    Field<Type>& iF = const_cast<Field<Type>&>(this->primitiveField());

    // Evaluate the intact patch field
    Field<Type> intactVals;
    intactPointPatchField_->evaluate(commsType);
    if (isA<valuePointPatchField<Type>>(intactPointPatchField_()))
    {
        intactVals =
            refCast<const valuePointPatchField<Type>>
            (
                intactPointPatchField_()
            );
    }
    else
    {
        intactVals = intactPointPatchField_->patchInternalField();
    }

    // Reset the internal field
    this->setInInternalField(iF, pIf);

    burstPointPatchField_->evaluate(commsType);

    // Calculate the mixed boundary
    (*this) =
        intactVals*burstBase_.pointIntact()
      + burstPointPatchField_->patchInternalField()
       *(1.0 - burstBase_.pointIntact());

    valuePointPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::burstPointPatchField<Type>::write(Ostream& os) const
{
    burstPointPatchField_->write(os);
    intactPointPatchField_->write(os);
    valuePointPatchField<Type>::write(os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::burstPointPatchField<Type>::operator=
(
    const burstPointPatchField<Type>& ptf
)
{
    burstPointPatchField_() = ptf.burstPatchField();
    intactPointPatchField_() = ptf.intactPatchField();
}


template<class Type>
void Foam::burstPointPatchField<Type>::operator=
(
    const pointPatchField<Type>& ptf
)
{
    burstPointPatchField_() = ptf;
    intactPointPatchField_() = ptf;
}


template<class Type>
void Foam::burstPointPatchField<Type>::operator=
(
    const Field<Type>& tf
)
{
    burstPointPatchField_() = tf;
    intactPointPatchField_() = tf;
}


template<class Type>
void Foam::burstPointPatchField<Type>::operator=
(
    const Type& t
)
{
    burstPointPatchField_() = t;
    intactPointPatchField_() = t;
}


template<class Type>
void Foam::burstPointPatchField<Type>::operator==
(
    const burstPointPatchField<Type>& ptf
)
{
    burstPointPatchField_() == ptf.intactPatchField();
    intactPointPatchField_() == ptf.intactPatchField();
}


template<class Type>
void Foam::burstPointPatchField<Type>::operator==
(
    const pointPatchField<Type>& ptf
)
{
    burstPointPatchField_() == ptf;
    intactPointPatchField_() == ptf;
}


template<class Type>
void Foam::burstPointPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    burstPointPatchField_() == tf;
    intactPointPatchField_() == tf;
}


template<class Type>
void Foam::burstPointPatchField<Type>::operator==
(
    const Type& t
)
{
    burstPointPatchField_() == t;
    intactPointPatchField_() == t;
}


// ************************************************************************* //

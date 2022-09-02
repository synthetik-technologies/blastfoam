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

#include "burstCyclicPointPatchField.H"
#include "calculatedPointPatchField.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::burstCyclicPointPatchField<Type>::burstCyclicPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    cyclicPointPatchField<Type>(p, iF),
    burstPointPatchFieldBase
    (
        dynamicCast<const facePointPatch>(this->patch()).patch()
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
Foam::burstCyclicPointPatchField<Type>::burstCyclicPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    cyclicPointPatchField<Type>(p, iF, dict),
    burstPointPatchFieldBase
    (
        dynamicCast<const facePointPatch>(this->patch()).patch()
    ),
    intactPointPatchField_()
{
    if (!isType<burstCyclicPointPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "patch " << this->patch().index() << " not burstCyclic type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }

    // Create a new patch dictionary and replace the type with the intactType
    dictionary intactDict(dict.parent(), dict);
    intactDict.set("type", dict.lookup<word>("intactType"));
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


template<class Type>
Foam::burstCyclicPointPatchField<Type>::burstCyclicPointPatchField
(
    const burstCyclicPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    cyclicPointPatchField<Type>(ptf, p, iF, mapper),
    burstPointPatchFieldBase
    (
        dynamicCast<const facePointPatch>(this->patch()).patch()
    ),
    intactPointPatchField_(ptf.intactPointPatchField_->clone(iF).ptr())
{
    if (!isType<burstCyclicPointPatch>(this->patch()))
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
Foam::burstCyclicPointPatchField<Type>::burstCyclicPointPatchField
(
    const burstCyclicPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    cyclicPointPatchField<Type>(ptf, iF),
    burstPointPatchFieldBase
    (
        dynamicCast<const facePointPatch>(this->patch()).patch()
    ),
    intactPointPatchField_(ptf.intactPointPatchField_->clone(iF).ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::burstCyclicPointPatchField<Type>::autoMap
(
    const pointPatchFieldMapper& m
)
{
    cyclicPointPatchField<Type>::autoMap(m);
    intactPointPatchField_->autoMap(m);
}


template<class Type>
void Foam::burstCyclicPointPatchField<Type>::rmap
(
    const pointPatchField<Type>& ptf,
    const labelList& addr
)
{
    cyclicPointPatchField<Type>::rmap(ptf, addr);

    const burstCyclicPointPatchField<Type>& bpf =
        refCast<const burstCyclicPointPatchField<Type>>(ptf);
    intactPointPatchField_->rmap(bpf.intactPointPatchField_(), addr);
}


template<class Type>
void Foam::burstCyclicPointPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    intactPointPatchField_->updateCoeffs();
    cyclicPointPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::burstCyclicPointPatchField<Type>::evaluate
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
        // Evaluate the cyclic patch
        intactPointPatchField_->evaluate(commsType);

        intactVals =
            refCast<const valuePointPatchField<Type>>
            (
                intactPointPatchField_()
            );
    }
    else
    {
        intactPointPatchField_->evaluate(commsType);
        intactVals = intactPointPatchField_->patchInternalField();
    }

    // Reset the internal field
    this->setInInternalField(iF, pIf);

    // Calculate the mixed boundary
    pIf =
        intactVals*burstBase_.pointIntact()
      + this->patchInternalField()*(1.0 - burstBase_.pointIntact());

    // Set the internal field
    this->setInInternalField(iF, pIf);

    pointPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::burstCyclicPointPatchField<Type>::write(Ostream& os) const
{
    {
        // Writing is a little weird since the intactPatchField has a different
        // type, but is in the same dictionary
        OStringStream oss;
        intactPointPatchField_->write(oss);
        dictionary dict(IStringStream(oss.str())());

        dict.changeKeyword("type", "intactType", false);
        dict.remove("patchType");
        forAllConstIter(IDLList<entry>, dict, iter)
        {
            iter().write(os);
        }
    }

    cyclicPointPatchField<Type>::write(os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::burstCyclicPointPatchField<Type>::operator=
(
    const burstCyclicPointPatchField<Type>& ptf
)
{
    intactPointPatchField_() = ptf.intactPatchField();
}


template<class Type>
void Foam::burstCyclicPointPatchField<Type>::operator=
(
    const pointPatchField<Type>& ptf
)
{
    intactPointPatchField_() = ptf;
}


template<class Type>
void Foam::burstCyclicPointPatchField<Type>::operator=
(
    const Field<Type>& tf
)
{
    intactPointPatchField_() = tf;
}


template<class Type>
void Foam::burstCyclicPointPatchField<Type>::operator=
(
    const Type& t
)
{
    intactPointPatchField_() = t;
}


template<class Type>
void Foam::burstCyclicPointPatchField<Type>::operator==
(
    const burstCyclicPointPatchField<Type>& ptf
)
{
    intactPointPatchField_() == ptf.intactPatchField();
}


template<class Type>
void Foam::burstCyclicPointPatchField<Type>::operator==
(
    const pointPatchField<Type>& ptf
)
{
    intactPointPatchField_() == ptf;
}


template<class Type>
void Foam::burstCyclicPointPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    intactPointPatchField_() == tf;
}


template<class Type>
void Foam::burstCyclicPointPatchField<Type>::operator==
(
    const Type& t
)
{
    intactPointPatchField_() == t;
}


// ************************************************************************* //

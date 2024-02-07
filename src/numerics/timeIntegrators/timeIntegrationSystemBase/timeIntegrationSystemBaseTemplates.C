/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

#include "timeIntegrationSystemBase.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FieldType>
void Foam::timeIntegrationSystemBase::storeOld
(
    FieldType& f,
    PtrList<FieldType>& fList,
    const bool conservative
)
{
    if (timeInt_->firstStep() && !timeInt_->restart())
    {
        f.storeOldTimes();
    }

    // Store fields if needed later
    label i = timeInt_->getOldIndex(step());
    if (i != -1)
    {
        if (fList.set(i))
        {
            fList[i] = f;
        }
        else
        {
            fList.set
            (
                i,
                new FieldType
                (
                    IOobject::groupName(f.name() + "_old", Foam::name(i)),
                    f
                )
            );
        }

        // Scale old field for mesh motion before storage
        // if conservative
        if (conservative)
        {
            tmp<scalarField> V0(timeInt_->V0());
            if (V0.valid())
            {
                fList[i].field() *= V0;
            }
        }
        fList[i].checkIn();
    }
}


template<class FieldType>
void Foam::timeIntegrationSystemBase::storeDelta
(
    const FieldType& f,
    PtrList<FieldType>& fList
)
{
    // Store fields if needed later
    label i = timeInt_->getDeltaIndex(step());
    if (i != -1)
    {
        if (fList.set(i))
        {
            fList[i] = f;
        }
        else
        {
            fList.set
            (
                i,
                new FieldType
                (
                    IOobject::groupName(f.name() + "_delta", Foam::name(i)),
                    f
                )
            );
            fList[i].checkIn();
        }
    }
}


template<class Type>
void Foam::timeIntegrationSystemBase::storeOld
(
    Type& f,
    List<Type>& fList,
    const bool conservative
)
{
    if (conservative)
    {
        f *= timeInt_->totalV0();
    }

    // Store fields if needed later
    const label i = timeInt_->getOldIndex(step());
    if (i >= 0)
    {
        fList[i] = f;
    }
}


template<class Type>
void Foam::timeIntegrationSystemBase::storeDelta
(
    const Type& f,
    List<Type>& fList
)
{
    // Store fields if needed later
    const label i = timeInt_->getDeltaIndex(step());
    if (i >= 0)
    {
        fList[i] = f;
    }
}


template<class FieldType>
void Foam::timeIntegrationSystemBase::blendOld
(
    FieldType& f,
    const PtrList<FieldType>& fList,
    const bool conservative
) const
{
    const scalarList& scales = a();
    const label curStep = timeInt_->step();

    tmp<scalarField> tV0;
    tmp<scalarField> tV;
    if (conservative)
    {
        tV0 = timeInt_->V0();
        tV = timeInt_->V();
    }

    // Scale current step by weight
    f *= scales[curStep];
    if (tV0.valid())
    {
        f.field() *= tV0();
    }

    forAll(scales, stepi)
    {
        if (curStep != stepi)
        {
            label i = timeInt_->getOldIndex(stepi);
            if (i != -1 && scales[stepi] != 0)
            {
                f += scales[stepi]*fList[i];
            }
        }
    }

    if (tV.valid())
    {
        f.field() /= tV();
    }
}


template<class Type>
void Foam::timeIntegrationSystemBase::blendOld
(
    Type& f,
    const List<Type>& fList,
    const bool conservative
) const
{
    const scalarList& scales = a();
    const label curStep = timeInt_->step();

    // Scale current step by weight
    f *= scales[curStep];
    if (conservative)
    {
        f *= timeInt_->totalV0();
    }

    forAll(scales, stepi)
    {
        if (curStep != stepi)
        {
            label i = timeInt_->getOldIndex(stepi);
            if (i != -1 && scales[stepi] != 0)
            {
                f += scales[stepi]*fList[i];
            }
        }
    }

    if (conservative)
    {
        f /= timeInt_->totalV();
    }
}


template<template<class> class ListType, class Type>
void Foam::timeIntegrationSystemBase::blendDelta
(
    Type& f,
    const ListType<Type>& fList
) const
{
    const scalarList& scales = this->b();
    const label curStep = timeInt_->step();

    // Scale current step by weight
    f *= scales[curStep];
    forAll(scales, stepi)
    {
        if (curStep != stepi)
        {
            label i = timeInt_->getDeltaIndex(stepi);
            if (i != -1 && scales[stepi] != 0)
            {
                f += scales[stepi]*fList[i];
            }
        }
    }
}


template<template<class> class ListType, class Type>
void Foam::timeIntegrationSystemBase::storeAndBlendOld
(
    Type& f,
    ListType<Type>& fList,
    const bool conservative
)
{
    storeOld(f, fList, conservative);
    blendOld(f, fList, conservative);
}


template<template<class> class ListType, class Type>
void Foam::timeIntegrationSystemBase::storeAndBlendDelta
(
    Type& f,
    ListType<Type>& fList
)
{
    storeDelta(f, fList);
    blendDelta(f, fList);
}


template<template<class> class ListType, class Type>
Foam::tmp<Type> Foam::timeIntegrationSystemBase::calcDelta
(
    const Type& f,
    const ListType<Type>& fList
) const
{
    tmp<Type> fN(new Type(f));
    const scalarList& scales = b();
    const label curStep = timeInt_->step();

    if (scales[curStep] == 0)
    {
        return fN*0.0;
    }

    // Remove old steps
    forAll(scales, stepi)
    {
        if (curStep != stepi)
        {
            label fi = timeInt_->getDeltaIndex(stepi);
            if (fi != -1 && scales[fi] != 0)
            {
                fN.ref() -= scales[fi]*fList[fi];
            }
        }
    }
    fN.ref() /= scales[curStep];
    return fN;
}


template<template<class> class ListType, class Type>
Foam::tmp<Type> Foam::timeIntegrationSystemBase::calcAndStoreDelta
(
    const Type& f,
    ListType<Type>& fList
)
{
    tmp<Type> fN(calcDelta(f, fList));
    storeDelta(fN(), fList);
    return fN;
}


// ************************************************************************* //

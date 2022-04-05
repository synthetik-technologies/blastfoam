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

#include "timeIntegrationSystem.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FieldType>
void Foam::timeIntegrationSystem::storeOld
(
    FieldType& f,
    PtrList<FieldType>& fList,
    const bool conservative
)
{
    if (step() == 1)
    {
        f.storeOldTimes();

        // Correct old field for mesh motion before storage
        if (meshPtr_->moving() && conservative)
        {
            f.ref() *= timeInt_->V0byV();
        }
    }

    // Store fields if needed later
    label i = oldIs_[step() - 1];
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
                FieldType::New
                (
                    f.name() + "_old_" + Foam::name(step() - 1),
                    f
                )
            );
        }
    }
}


template<class FieldType>
void Foam::timeIntegrationSystem::storeOld(FieldType& f, const bool conservative)
{
    if (!timeInt_.valid())
    {
        return;
    }
    storeOld(f, timeInt_->oldFields(f), conservative);
}


template<class FieldType>
void Foam::timeIntegrationSystem::storeDelta
(
    const FieldType& f,
    PtrList<FieldType>& fList
)
{
    // Store fields if needed later
    label i = deltaIs_[step() - 1];
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
                FieldType::New
                (
                    f.name() + "_delta_" + Foam::name(step() - 1),
                    f
                )
            );
        }
    }
}


template<class FieldType>
void Foam::timeIntegrationSystem::storeDelta(const FieldType& f)
{
    if (!timeInt_.valid())
    {
        return;
    }
    storeDelta(f, timeInt_->deltaFields(f));
}


template<class Type>
void Foam::timeIntegrationSystem::storeOld
(
    Type& f,
    List<Type>& fList,
    const bool
)
{
    // Store fields if needed later
    if (oldIs_[step() - 1] != -1)
    {
        fList[oldIs_[step() - 1]] = f;
    }
}


template<class Type>
void Foam::timeIntegrationSystem::storeDelta
(
    const Type& f,
    List<Type>& fList
)
{
    // Store fields if needed later
    if (deltaIs_[step() - 1] != -1)
    {
        fList[deltaIs_[step() - 1]] = f;
    }
}


template<template<class> class ListType, class Type>
void Foam::timeIntegrationSystem::blendOld
(
    Type& f,
    const ListType<Type>& fList
) const
{
    blendSteps(oldIs_, f, fList, a());
}


template<class FieldType>
void Foam::timeIntegrationSystem::blendOld(FieldType& f) const
{
    if (!timeInt_.valid())
    {
        return;
    }
    blendSteps(oldIs_, f, timeInt_->oldFields(f), a());
}


template<template<class> class ListType, class Type>
void Foam::timeIntegrationSystem::blendDelta
(
    Type& f,
    const ListType<Type>& fList
) const
{
    blendSteps(deltaIs_, f, fList, b());
}


template<class FieldType>
void Foam::timeIntegrationSystem::blendDelta(FieldType& f) const
{
    if (!timeInt_.valid())
    {
        return;
    }
    blendSteps(oldIs_, f, timeInt_->deltaFields(f), b());
}


template<template<class> class ListType, class Type>
void Foam::timeIntegrationSystem::storeAndBlendOld
(
    Type& f,
    ListType<Type>& fList,
    const bool conservative
)
{
    storeOld(f, fList, conservative);
    blendSteps(oldIs_, f, fList, a());
}


template<class FieldType>
void Foam::timeIntegrationSystem::storeAndBlendOld
(
    FieldType& f,
    const bool conservative
)
{
    if (!timeInt_.valid())
    {
        return;
    }
    storeAndBlendOld(f, timeInt_->oldFields(f), conservative);
}


template<template<class> class ListType, class Type>
void Foam::timeIntegrationSystem::storeAndBlendDelta
(
    Type& f,
    ListType<Type>& fList
)
{
    storeDelta(f, fList);
    blendSteps(deltaIs_, f, fList, b());
}

template<class FieldType>
void Foam::timeIntegrationSystem::storeAndBlendDelta(FieldType& f)
{
    if (!timeInt_.valid())
    {
        return;
    }
    storeAndBlendDelta(f, timeInt_->deltaFields(f));
}

template<template<class> class ListType, class Type>
Foam::tmp<Type> Foam::timeIntegrationSystem::calcDelta
(
    const Type& f,
    const ListType<Type>& fList
) const
{
    tmp<Type> fN(new Type(f));
    const scalarList& scales(b());
    const labelList& indices(deltaIs_);

    if (scales[step() - 1] == 0)
    {
        return fN*0.0;
    }

    // Remove old steps
    for (label i = 0; i < step() - 1; i++)
    {
        label fi = indices[i];
        if (fi != -1 && scales[fi] != 0)
        {
            fN.ref() -= scales[fi]*fList[fi];
        }
    }
    fN.ref() /= scales[step() - 1];
    return fN;
}


template<class FieldType>
Foam::tmp<FieldType> Foam::timeIntegrationSystem::calcDelta(const FieldType& f) const
{
    if (!timeInt_.valid())
    {
        return FieldType::New(f.name(), f);
    }
    return calcDelta(f, timeInt_->deltaFields(f));
}


template<template<class> class ListType, class Type>
Foam::tmp<Type> Foam::timeIntegrationSystem::calcAndStoreDelta
(
    const Type& f,
    ListType<Type>& fList
)
{
    tmp<Type> fN(calcDelta(f, fList));
    storeDelta(fN(), fList);
    return fN;
}


template<class FieldType>
Foam::tmp<FieldType> Foam::timeIntegrationSystem::calcAndStoreDelta(const FieldType& f)
{
    if (!timeInt_.valid())
    {
        return tmp<FieldType>(FieldType::New(f.name(), f));
    }
   return calcAndStoreDelta(f, timeInt_->deltaFields(f));
}


template<template<class> class ListType, class Type>
void Foam::timeIntegrationSystem::blendSteps
(
    const labelList& indices,
    Type& f,
    const ListType<Type>& fList,
    const scalarList& scales
) const
{
    // Scale current step by weight
    f *= scales[step() - 1];
    for (label i = 0; i < step() - 1; i++)
    {
        label fi = indices[i];
        if (fi != -1 && scales[fi] != 0)
        {
            f += scales[fi]*fList[fi];
        }
    }
}


template<class FieldType>
void Foam::timeIntegrationSystem::addOldField(const FieldType& f)
{
    if (!timeInt_.valid())
    {
        return;
    }
    timeInt_->addOldField(f);
}


template<class FieldType>
void Foam::timeIntegrationSystem::addDeltaField(const FieldType& f)
{
    if (!timeInt_.valid())
    {
        return;
    }
    timeInt_->addDeltaField(f);
}


template<class FieldType>
void Foam::timeIntegrationSystem::addOldDeltaField(const FieldType& f)
{
    if (!timeInt_.valid())
    {
        return;
    }
    timeInt_->addOldField(f);
    timeInt_->addDeltaField(f);
}

// ************************************************************************* //

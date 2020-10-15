/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2020
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

#include "integrationSystem.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class fieldType>
void Foam::integrationSystem::storeOld
(
    fieldType& f,
    PtrList<fieldType>& fList,
    const bool moving
) const
{
    // Correct old field for mesh motion before storage
    if (f.mesh().moving() && step() == 1 && moving)
    {
        f.ref() *= f.mesh().V0()/f.mesh().V();
    }

    // Store fields if needed later
    if (oldIs_[step() - 1] != -1)
    {
        fList.set
        (
            oldIs_[step() - 1],
            new fieldType(f.name() + Foam::name(step() - 1), f)
        );
    }
}


template<class fieldType>
void Foam::integrationSystem::storeDelta
(
    const fieldType& f,
    PtrList<fieldType>& fList
) const
{
    // Store fields if needed later
    if (deltaIs_[step() - 1] != -1)
    {
        fList.set
        (
            deltaIs_[step() - 1],
            new fieldType(f.name() + Foam::name(step() - 1), f)
        );
    }
}


template<class Type>
void Foam::integrationSystem::storeOld
(
    Type& f,
    List<Type>& fList,
    const bool
) const
{
    // Store fields if needed later
    if (oldIs_[step() - 1] != -1)
    {
        fList[oldIs_[step() - 1]] = f;
    }
}


template<class Type>
void Foam::integrationSystem::storeDelta
(
    const Type& f,
    List<Type>& fList
) const
{
    // Store fields if needed later
    if (deltaIs_[step() - 1] != -1)
    {
        fList[deltaIs_[step() - 1]] = f;
    }
}


template<template<class> class ListType, class Type>
void Foam::integrationSystem::blendOld
(
    Type& f,
    const ListType<Type>& fList
) const
{
    blendSteps(oldIs_, f, fList, a());
}


template<template<class> class ListType, class Type>
void Foam::integrationSystem::blendDelta
(
    Type& f,
    const ListType<Type>& fList
) const
{
    blendSteps(deltaIs_, f, fList, b());
}


template<template<class> class ListType, class Type>
void Foam::integrationSystem::blendSteps
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


template<template<class> class ListType, class Type>
void Foam::integrationSystem::storeAndBlendOld
(
    Type& f,
    ListType<Type>& fList,
    const bool moving
) const
{
    storeOld(f, fList, moving);
    blendSteps(oldIs_, f, fList, a());
}



template<template<class> class ListType, class Type>
void Foam::integrationSystem::storeAndBlendDelta
(
    Type& f,
    ListType<Type>& fList
) const
{
    storeDelta(f, fList);
    blendSteps(deltaIs_, f, fList, b());
}


template<class fieldType>
void Foam::integrationSystem::clearOld(PtrList<fieldType>& fList) const
{
    fList.clear();
    fList.resize(nOld_);
}


template<class fieldType>
void Foam::integrationSystem::clearDelta(PtrList<fieldType>& fList) const
{
    fList.clear();
    fList.resize(nDelta_);
}
// ************************************************************************* //

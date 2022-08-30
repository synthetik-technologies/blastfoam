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

#include "timeIntegrator.H"

// * * * * * * * * * * * * * * Protected Functions * * * * * * * * * * * * * //

template<class FieldType>
const Foam::PtrList<FieldType>& Foam::timeIntegrator::lookupOld
(
    HashPtrTable<PtrList<FieldType>>& table,
    const FieldType& f
) const
{
    if (!table.found(f.name()))
    {
        insertOldList(table, f);
    }
    return *table[f.name()];
}


template<class FieldType>
Foam::PtrList<FieldType>& Foam::timeIntegrator::lookupOld
(
    HashPtrTable<PtrList<FieldType>>& table,
    const FieldType& f
)
{
    if (!table.found(f.name()))
    {
        insertOldList(table, f);
    }
    return *table[f.name()];
}


template<class FieldType>
const Foam::PtrList<FieldType>& Foam::timeIntegrator::lookupDelta
(
    HashPtrTable<PtrList<FieldType>>& table,
    const FieldType& f
) const
{
    if (!table.found(f.name()))
    {
        insertDeltaList(table, f);
    }
    return *table[f.name()];
}


template<class FieldType>
Foam::PtrList<FieldType>& Foam::timeIntegrator::lookupDelta
(
    HashPtrTable<PtrList<FieldType>>& table,
    const FieldType& f
)
{
    if (!table.found(f.name()))
    {
        insertDeltaList(table, f);
    }
    return *table[f.name()];
}


template<class FieldType>
void Foam::timeIntegrator::insertOldList
(
    HashPtrTable<PtrList<FieldType>>& table,
    const FieldType& f
) const
{
    table.insert(f.name(), new PtrList<FieldType>(nOld_));
}

//- Insert a list of delta fields into the given hash table
template<class FieldType>
void Foam::timeIntegrator::insertDeltaList
(
    HashPtrTable<PtrList<FieldType>>& table,
    const FieldType& f
) const
{
    table.insert(f.name(), new PtrList<FieldType>(nDelta_));
}


template<class FieldType>
void Foam::timeIntegrator::clearOldFields
(
    HashPtrTable<PtrList<FieldType>>& table
)
{
    forAllIter
    (
        typename HashPtrTable<PtrList<FieldType>>,
        table,
        iter
    )
    {
        savedOldFields_.insert(iter.key());
        label size = iter()->size();
        iter()->clear();
        iter()->resize(size);
    }
}


template<class FieldType>
void Foam::timeIntegrator::clearDeltaFields
(
    HashPtrTable<PtrList<FieldType>>& table
)
{
    forAllIter
    (
        typename HashPtrTable<PtrList<FieldType>>,
        table,
        iter
    )
    {
        label size = iter()->size();
        iter()->clear();
        iter()->resize(size);
    }
}


template<class FieldType>
void Foam::timeIntegrator::resetFields()
{
    forAllConstIter(wordHashSet, savedOldFields_, iter)
    {
        if (mesh_.foundObject<FieldType>(iter.key()))
        {
            FieldType& f = mesh_.lookupObjectRef<FieldType>(iter.key());
            f == f.oldTime();
        }
    }
}

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
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

#include "fvTimeIntegrator.H"

// * * * * * * * * * * * * * * Protected Functions * * * * * * * * * * * * * //

template<class FieldType>
const Foam::PtrList<FieldType>& Foam::fvTimeIntegrator::lookupOld
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
Foam::PtrList<FieldType>& Foam::fvTimeIntegrator::lookupOldRef
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
const Foam::PtrList<FieldType>& Foam::fvTimeIntegrator::lookupDelta
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
Foam::PtrList<FieldType>& Foam::fvTimeIntegrator::lookupDeltaRef
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
void Foam::fvTimeIntegrator::insertOldList
(
    HashPtrTable<PtrList<FieldType>>& table,
    const FieldType& f
) const
{
    table.insert(f.name(), new PtrList<FieldType>(nOld_));
}

//- Insert a list of delta fields into the given hash table
template<class FieldType>
void Foam::fvTimeIntegrator::insertDeltaList
(
    HashPtrTable<PtrList<FieldType>>& table,
    const FieldType& f
) const
{
    table.insert(f.name(), new PtrList<FieldType>(nDelta_));
}


template<class FieldType>
void Foam::fvTimeIntegrator::clearOldFields
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
        iter()->clear();
        iter()->resize(nOld_);
    }
}


template<class FieldType>
void Foam::fvTimeIntegrator::clearDeltaFields
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
        iter()->clear();
        iter()->resize(nDelta_);
    }
}


template<class FieldType>
void Foam::fvTimeIntegrator::resetFields()
{
    forAllConstIter(wordHashSet, savedOldFields_, iter)
    {
        if (mesh_.foundObject<FieldType>(iter.key()))
        {
            DebugInfo<< "Resetting " << iter.key() << endl;
            FieldType& f = mesh_.lookupObjectRef<FieldType>(iter.key());
            f == f.oldTime();
        }
    }
}


template<class FieldType>
void Foam::fvTimeIntegrator::shuffleFields
(
    HashPtrTable<PtrList<FieldType>>& table,
    const labelList& oldToNew,
    const word& type
)
{
    forAllIter
    (
        typename HashPtrTable<PtrList<FieldType>>,
        table,
        iter
    )
    {
        // Make sure everything is the correct size
        PtrList<FieldType>& fields = *iter();

        // Check out all fields to make sure we can rename
        forAll(fields, i)
        {
            if (fields.set(i))
            {
                fields[i].checkOut();
            }
        }
        fields.setSize(oldToNew.size());


        // Temporarily transfer fields to a new list so shuffling does not overwrite data
        PtrList<FieldType> tmp(fields.size());

        // Reorder based on the mapping
        forAll(oldToNew, i)
        {
            const label newi = oldToNew[i];
            if (newi >= 0)
            {
                tmp.set(newi, fields.set(i, nullptr).ptr());
                if (tmp.set(newi))
                {
                    tmp[newi].rename
                    (
                        IOobject::groupName
                        (
                            iter.key() + "_" + type,
                            Foam::name(newi)
                        )
                    );
                    tmp[newi].checkIn();
                }
            }
        }
        fields.transfer(tmp);
    }
}


// ************************************************************************* //

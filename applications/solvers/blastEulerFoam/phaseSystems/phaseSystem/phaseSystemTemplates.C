/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class modelType>
void Foam::phaseSystem::createSubModels
(
    const dictTable& modelDicts,
    HashTable
    <
        autoPtr<modelType>,
        phasePairKey,
        phasePairKey::hash
    >& models
)
{
    forAllConstIter(dictTable, modelDicts, iter)
    {
        const phasePairKey& key = iter.key();

        models.insert
        (
            key,
            modelType::New
            (
               *iter,
                phasePairs_[key]
            )
        );
    }
}


template<class modelType>
void Foam::phaseSystem::generatePairsAndSubModels
(
    const word& modelName,
    HashTable
    <
        autoPtr<modelType>,
        phasePairKey,
        phasePairKey::hash
    >& models
)
{
    dictTable modelDicts(lookup(modelName));

    generatePairs(modelDicts);

    createSubModels(modelDicts, models);
}


template<class modelType>
void Foam::phaseSystem::generatePairsAndSubModels
(
    const word& modelName,
    HashTable
    <
        autoPtr<modelType>,
        phasePairKey,
        phasePairKey::hash
    >& models,
    const bool correctFixedFluxBCs
)
{
    typedef
        HashTable<autoPtr<modelType>, phasePairKey, phasePairKey::hash>
        modelTypeTable;

    modelTypeTable tempModels;
    generatePairsAndSubModels(modelName, tempModels);

    autoPtr<modelType> noModel(nullptr);

    forAllConstIter(typename modelTypeTable, tempModels, iter)
    {
        if (!iter().valid())
        {
            continue;
        }

        const phasePairKey key(iter.key().first(), iter.key().second());


        models.insert
        (
            key,
            autoPtr<modelType>
            (
                new modelType
                (
                    tempModels.found(key) ? tempModels[key] : noModel
                )
            )
        );

        if (!phasePairs_.found(key))
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new phasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class modelType>
bool Foam::phaseSystem::foundSubModel(const phasePair& key) const
{
    const word name(IOobject::groupName(modelType::typeName, key.name()));

    if (key.ordered())
    {
        if (mesh().foundObject<modelType>(name))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        if
        (
            mesh().foundObject<modelType>(name)
         ||
            mesh().foundObject<modelType>
            (
                IOobject::groupName(modelType::typeName, key.otherName())
            )
        )
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}


template<class modelType>
const modelType& Foam::phaseSystem::lookupSubModel(const phasePair& key) const
{
    const word name(IOobject::groupName(modelType::typeName, key.name()));

    if (key.ordered() || mesh().foundObject<modelType>(name))
    {
        return mesh().lookupObject<modelType>(name);
    }
    else
    {
        return
            mesh().lookupObject<modelType>
            (
                IOobject::groupName(modelType::typeName, key.otherName())
            );
    }
}


template<class modelType>
bool Foam::phaseSystem::foundSubModel
(
    const phaseModel& dispersed,
    const phaseModel& continuous
) const
{
    return foundSubModel<modelType>(orderedPhasePair(dispersed, continuous));
}


template<class modelType>
const modelType& Foam::phaseSystem::lookupSubModel
(
    const phaseModel& dispersed,
    const phaseModel& continuous
) const
{
    return lookupSubModel<modelType>(orderedPhasePair(dispersed, continuous));
}


// ************************************************************************* //

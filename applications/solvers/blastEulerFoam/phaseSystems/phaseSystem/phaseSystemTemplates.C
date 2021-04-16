/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

#include "BlendedInterfacialModel.H"

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
    >& models,
    const bool required
)
{
    if (!found(modelName))
    {
        return;
    }

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
        autoPtr<BlendedInterfacialModel<modelType>>,
        phasePairKey,
        phasePairKey::hash
    >& models,
    const bool required,
    const bool correctFixedFluxBCs
)
{
    typedef
        HashTable<autoPtr<modelType>, phasePairKey, phasePairKey::hash>
        modelTypeTable;

    modelTypeTable tempModels;
    generatePairsAndSubModels(modelName, tempModels, false);

    const blendingMethod& blending
    (
        blendingMethods_.found(modelName)
      ? blendingMethods_[modelName]
      : blendingMethods_.found(member(modelName))
      ? blendingMethods_[member(modelName)]
      : blendingMethods_["default"]
    );

    autoPtr<modelType> noModel(nullptr);

    forAllConstIter(typename modelTypeTable, tempModels, iter)
    {
        if (!iter().valid())
        {
            continue;
        }

        const phasePairKey key(iter.key().first(), iter.key().second());
        const phasePairKey key1In2(key.first(), key.second(), true);
        const phasePairKey key2In1(key.second(), key.first(), true);

        models.insert
        (
            key,
            autoPtr<BlendedInterfacialModel<modelType>>
            (
                new BlendedInterfacialModel<modelType>
                (
                    phaseModels_[key.first()],
                    phaseModels_[key.second()],
                    blending,
                    tempModels.found(key    ) ? tempModels[key    ] : noModel,
                    tempModels.found(key1In2) ? tempModels[key1In2] : noModel,
                    tempModels.found(key2In1) ? tempModels[key2In1] : noModel,
                    correctFixedFluxBCs
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


template<class modelType>
void Foam::phaseSystem::generatePairsAndSubModels
(
    const word& modelName,
    HashTable
    <
        Pair<autoPtr<modelType>>,
        phasePairKey,
        phasePairKey::hash
    >& models,
    const bool require,
    const bool correctFixedFluxBCs
)
{
    typedef
        HashTable<autoPtr<modelType>, phasePairKey, phasePairKey::hash>
        modelTypeTable;

    forAll(phaseModels_, phasei)
    {
        const phaseModel& phase = phaseModels_[phasei];

        modelTypeTable tempModels;
        generatePairsAndSubModels
        (
            IOobject::groupName(modelName, phase.name()),
            tempModels,
            require,
            correctFixedFluxBCs
        );

        forAllIter(typename modelTypeTable, tempModels, tempModelIter)
        {
            const phasePairKey& key(tempModelIter.key());

            if (!models.found(key))
            {
                models.insert
                (
                    key,
                    Pair<autoPtr<modelType>>()
                );
            }

            const phasePair& pair = phasePairs_[key];

            if (!pair.contains(phase))
            {
                FatalErrorInFunction
                    << "A two-sided " << modelType::typeName << " was "
                    << "specified for the " << phase.name() << " side of the "
                    << pair << " pair, but that phase is not part of that pair."
                    << exit(FatalError);
            }

            models[key][pair.index(phase)] = tempModelIter().ptr();
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
    const phaseModel& continuous,
    const bool ordered
) const
{
    if (ordered)
    {
        return foundSubModel<modelType>(orderedPhasePair(dispersed, continuous));
    }
    return foundSubModel<modelType>(phasePair(dispersed, continuous));
}


template<class modelType>
const modelType& Foam::phaseSystem::lookupSubModel
(
    const phaseModel& dispersed,
    const phaseModel& continuous,
    const bool ordered
) const
{
    if (ordered)
    {
        return lookupSubModel<modelType>(orderedPhasePair(dispersed, continuous));
    }
    return lookupSubModel<modelType>(phasePair(dispersed, continuous));
}


template<class modelType>
bool Foam::phaseSystem::foundBlendedSubModel(const phasePair& key) const
{
    if
    (
        mesh().foundObject<BlendedInterfacialModel<modelType>>
        (
            IOobject::groupName
            (
                BlendedInterfacialModel<modelType>::typeName,
                key.name()
            )
        )
     || mesh().foundObject<BlendedInterfacialModel<modelType>>
        (
            IOobject::groupName
            (
                BlendedInterfacialModel<modelType>::typeName,
                key.otherName()
            )
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


template<class modelType>
const Foam::BlendedInterfacialModel<modelType>&
Foam::phaseSystem::lookupBlendedSubModel(const phasePair& key) const
{
    const word name
    (
        IOobject::groupName
        (
            BlendedInterfacialModel<modelType>::typeName,
            key.name()
        )
    );

    if (mesh().foundObject<BlendedInterfacialModel<modelType>>(name))
    {
        return mesh().lookupObject<BlendedInterfacialModel<modelType>>(name);
    }
    else
    {
        return
            mesh().lookupObject<BlendedInterfacialModel<modelType>>
            (
                IOobject::groupName
                (
                    BlendedInterfacialModel<modelType>::typeName,
                    key.otherName()
                )
            );
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Synthetik Applied Technologies:  |  Store field maxes over all times
07-09-2021 Synthetik Applied Technologies:  |  Added minimum option
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "fieldMinMax.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldMin, 0);
    defineTypeNameAndDebug(fieldMax, 0);
    defineTypeNameAndDebug(fieldMinMax, 0);
    addToRunTimeSelectionTable(functionObject, fieldMin, dictionary);
    addToRunTimeSelectionTable(functionObject, fieldMax, dictionary);
    addToRunTimeSelectionTable(functionObject, fieldMinMax, dictionary);
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldMinMax::modeType,
    3
>::names[] = {"component", "magnitude", "componentMag"};

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldMinMax::minMaxType,
    2
>::names[] = {"min", "max"};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldMinMax::modeType,
    3
> Foam::functionObjects::fieldMinMax::modeTypeNames_;

const Foam::NamedEnum
<
    Foam::functionObjects::fieldMinMax::minMaxType,
    2
> Foam::functionObjects::fieldMinMax::minMaxTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMinMax::fieldMinMax
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    restartOnRestart_(dict.lookupOrDefault("restartOnRestart", false)),
    mode_(modeType::cmpt),
    minMax_(minMaxTypeNames_[dict.lookupOrDefault<word>("minMax", "max")]),
    minMaxName_(minMax_ == minMaxType::min ? "Min" : "Max"),
    fieldNames_(dict.lookup("fields"))
{
    read(dict);
}

Foam::functionObjects::fieldMinMax::fieldMinMax
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const minMaxType mm
)
:
    fvMeshFunctionObject(name, runTime, dict),
    restartOnRestart_(dict.lookupOrDefault("restartOnRestart", false)),
    mode_(modeType::cmpt),
    minMax_(mm),
    minMaxName_(minMax_ == minMaxType::min ? "Min" : "Max"),
    fieldNames_(dict.lookup("fields"))
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMinMax::~fieldMinMax()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldMinMax::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Log << type() << " " << name() << ":" << nl;

    dict.readIfPresent("restartOnRestart", restartOnRestart_);
    mode_ = modeTypeNames_[dict.lookupOrDefault<word>("mode", "component")];    minMaxName_ = minMax_ == minMaxType::min ? "Min" : "Max";

    Log << endl;

    return true;
}


void Foam::functionObjects::fieldMinMax::topoChange
(
    const polyTopoChangeMap& map
)
{
    forAll(fieldNames_, fieldi)
    {
        bool found = false;
        #define MapFields(Type, Patch, Mesh)                \
        found =                                             \
            found                                           \
         || mapField<GeometricField<Type, Patch, Mesh>>     \
            (                                               \
                fieldNames_[fieldi],                        \
                map                                         \
            );

        FOR_ALL_FIELD_TYPES(MapFields, fvPatchField, volMesh);
        FOR_ALL_FIELD_TYPES(MapFields, fvsPatchField, surfaceMesh);

        #undef MapFields

        if (!found)
        {
            cannotFindObject(fieldNames_[fieldi]);
        }
    }

    forAll(fieldNames_, fieldi)
    {
        bool found = false;
        #define StoreOldFields(Type, Patch, Mesh)           \
        found =                                             \
            found                                           \
            || storeOld<GeometricField<Type, Patch, Mesh>>     \
            (                                               \
                fieldNames_[fieldi]                         \
            );

        FOR_ALL_FIELD_TYPES(StoreOldFields, fvPatchField, volMesh);
        FOR_ALL_FIELD_TYPES(StoreOldFields, fvsPatchField, surfaceMesh);

        #undef StoreOldFields

        if (!found)
        {
            cannotFindObject(fieldNames_[fieldi]);
        }
    }

}


bool Foam::functionObjects::fieldMinMax::execute()
{
    forAll(fieldNames_, fieldi)
    {
        bool found = false;
        #define UpdateFields(Type, Patch, Mesh)             \
        found =                                             \
            found                                           \
         || update<GeometricField<Type, Patch, Mesh>>       \
            (                                               \
                fieldNames_[fieldi]                         \
            );

        FOR_ALL_FIELD_TYPES(UpdateFields, fvPatchField, volMesh);
        FOR_ALL_FIELD_TYPES(UpdateFields, fvsPatchField, surfaceMesh);

        #undef UpdateFields

        if (!found)
        {
            cannotFindObject(fieldNames_[fieldi]);
        }
    }

    if (mesh_.dynamic())
    {
        forAll(fieldNames_, fieldi)
        {
            bool found = false;
            #define StoreOldFields(Type, Patch, Mesh)           \
            found =                                             \
                found                                           \
             || storeOld<GeometricField<Type, Patch, Mesh>>     \
                (                                               \
                    fieldNames_[fieldi]                         \
                );

            FOR_ALL_FIELD_TYPES(StoreOldFields, fvPatchField, volMesh);
            FOR_ALL_FIELD_TYPES(StoreOldFields, fvsPatchField, surfaceMesh);

            #undef StoreOldFields

            if (!found)
            {
                cannotFindObject(fieldNames_[fieldi]);
            }
        }
    }
    return true;
}


bool Foam::functionObjects::fieldMinMax::write()
{
    if (obr_.time().timeIndex() == obr_.time().startTimeIndex())
    {
        return true;
    }

    bool good = true;
    forAll(fieldNames_, fieldi)
    {
        good = good && writeObject(computedName(fieldNames_[fieldi]));
    }
    return good;
}


// ************************************************************************* //

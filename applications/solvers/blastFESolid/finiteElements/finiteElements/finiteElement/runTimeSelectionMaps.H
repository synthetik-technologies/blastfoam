/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Description
    Macros to ease declaration of run-time selection maps.

    declareRunTimeSelectionMap is used to create a run-time selection map
    for a base-class which holds constructor pointers on the map.

    declareRunTimeNewSelectionMap is used to create a run-time selection
    map for a derived-class which holds "New" pointers on the map.

\*---------------------------------------------------------------------------*/

#ifndef runTimeSelectionMaps_H
#define runTimeSelectionMaps_H

#include "autoPtr.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Declare a run-time selection
#define declareRunTimeSelectionMap(autoPtr,baseType,argNames,argList,parList)   \
                                                                               \
    /* Construct from argList function pointer type */                         \
    typedef autoPtr<baseType> (*argNames##ConstructorPtr)argList;              \
                                                                               \
    /* Construct from argList function map type */                             \
    typedef Map<argNames##ConstructorPtr> argNames##ConstructorMap;            \
                                                                               \
    /* Construct from argList function pointer map pointer */                  \
    static argNames##ConstructorMap* argNames##ConstructorMapPtr_;             \
                                                                               \
    /* Map constructor called from the map add function */                     \
    static void construct##argNames##ConstructorMaps();                        \
                                                                               \
    /* Map destructor called from the map add function destructor */           \
    static void destroy##argNames##ConstructorMaps();                          \
                                                                               \
    /* Class to add constructor from argList to map */                         \
    template<class baseType##Type>                                             \
    class add##argNames##ConstructorToMap                                      \
    {                                                                          \
    public:                                                                    \
                                                                               \
        static autoPtr<baseType> New argList                                   \
        {                                                                      \
            return autoPtr<baseType>(new baseType##Type parList);              \
        }                                                                      \
                                                                               \
        add##argNames##ConstructorToMap(const label lookup)                    \
        {                                                                      \
            construct##argNames##ConstructorMaps();                            \
            if (!argNames##ConstructorMapPtr_->insert(lookup, New))            \
            {                                                                  \
                std::cerr<< "Duplicate entry " << lookup                       \
                    << " in runtime selection map " << #baseType               \
                    << std::endl;                                              \
                error::safePrintStack(std::cerr);                              \
            }                                                                  \
        }                                                                      \
                                                                               \
        ~add##argNames##ConstructorToMap()                                     \
        {                                                                      \
            destroy##argNames##ConstructorMaps();                              \
        }                                                                      \
    };


//- Declare a run-time selection for derived classes
#define declareRunTimeNewSelectionMap(                                         \
    autoPtr,baseType,argNames,argList,parList)                                 \
                                                                               \
    /* Construct from argList function pointer type */                         \
    typedef autoPtr<baseType> (*argNames##ConstructorPtr)argList;              \
                                                                               \
    /* Construct from argList function map type */                             \
    typedef Map<argNames##ConstructorPtr> argNames##ConstructorMap;            \
                                                                               \
    /* Construct from argList function pointer map pointer */                  \
    static argNames##ConstructorMap* argNames##ConstructorMapPtr_;             \
                                                                               \
    /* Map constructor called from the map add function */                     \
    static void construct##argNames##ConstructorMaps();                        \
                                                                               \
    /* Map destructor called from the map add function destructor */           \
    static void destroy##argNames##ConstructorMaps();                          \
                                                                               \
    /* Class to add constructor from argList to map */                         \
    template<class baseType##Type>                                             \
    class add##argNames##ConstructorToMap                                      \
    {                                                                          \
    public:                                                                    \
                                                                               \
        static autoPtr<baseType> New##baseType argList                         \
        {                                                                      \
            return autoPtr<baseType>(baseType##Type::New parList.ptr());       \
        }                                                                      \
                                                                               \
        add##argNames##ConstructorToMap(const label lookup)                    \
        {                                                                      \
            construct##argNames##ConstructorMaps();                            \
            if                                                                 \
            (                                                                  \
               !argNames##ConstructorMapPtr_->insert                           \
                (                                                              \
                    lookup,                                                    \
                    New##baseType                                              \
                )                                                              \
            )                                                                  \
            {                                                                  \
                std::cerr<< "Duplicate entry " << lookup                       \
                    << " in runtime selection map " << #baseType               \
                    << std::endl;                                              \
                error::safePrintStack(std::cerr);                              \
            }                                                                  \
        }                                                                      \
                                                                               \
        ~add##argNames##ConstructorToMap()                                     \
        {                                                                      \
            destroy##argNames##ConstructorMaps();                              \
        }                                                                      \
    };                                                                         \



// Constructor aid
#define defineRunTimeSelectionMapConstructor(baseType,argNames)                \
                                                                               \
    /* Map constructor called from the map add function */                     \
    void baseType::construct##argNames##ConstructorMaps()                      \
    {                                                                          \
        static bool constructed = false;                                       \
        if (!constructed)                                                      \
        {                                                                      \
            constructed = true;                                                \
            baseType::argNames##ConstructorMapPtr_                             \
                = new baseType::argNames##ConstructorMap;                      \
        }                                                                      \
    }


// Destructor aid
#define defineRunTimeSelectionMapDestructor(baseType,argNames)                 \
                                                                               \
    /* Map destructor called from the map add function destructor */           \
    void baseType::destroy##argNames##ConstructorMaps()                        \
    {                                                                          \
        if (baseType::argNames##ConstructorMapPtr_)                            \
        {                                                                      \
            delete baseType::argNames##ConstructorMapPtr_;                     \
            baseType::argNames##ConstructorMapPtr_ = nullptr;                  \
        }                                                                      \
    }


// Create pointer to hash-map of functions
#define defineRunTimeSelectionMapPtr(baseType,argNames)                        \
                                                                               \
    /* Define the constructor function map */                                  \
    baseType::argNames##ConstructorMap*                                        \
        baseType::argNames##ConstructorMapPtr_ = nullptr



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Define run-time selection map
#define defineRunTimeSelectionMap(baseType,argNames)                           \
                                                                               \
    defineRunTimeSelectionMapPtr(baseType,argNames);                           \
    defineRunTimeSelectionMapConstructor(baseType,argNames);                   \
    defineRunTimeSelectionMapDestructor(baseType,argNames)


//- Define run-time selection map for template classes
//  use when baseType doesn't need a template argument (eg, is a typedef)
#define defineTemplateRunTimeSelectionMap(baseType,argNames)                   \
                                                                               \
    template<>                                                                 \
    defineRunTimeSelectionMapPtr(baseType,argNames);                           \
    template<>                                                                 \
    defineRunTimeSelectionMapConstructor(baseType,argNames);                   \
    template<>                                                                 \
    defineRunTimeSelectionMapDestructor(baseType,argNames)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Constructor aid: use when baseType requires the Targ template argument
#define defineTemplatedRunTimeSelectionMapConstructor(baseType,argNames,Targ)  \
                                                                               \
    /* Map constructor called from the map add function */                     \
    void baseType<Targ>::construct##argNames##ConstructorMaps()                \
    {                                                                          \
        static bool constructed = false;                                       \
        if (!constructed)                                                      \
        {                                                                      \
            constructed = true;                                                \
            baseType<Targ>::argNames##ConstructorMapPtr_                       \
                = new baseType<Targ>::argNames##ConstructorMap;                \
        }                                                                      \
    }


// Destructor aid: use when baseType requires the Targ template argument
#define defineTemplatedRunTimeSelectionMapDestructor(baseType,argNames,Targ)   \
                                                                               \
    /* Map destructor called from the map add function destructor */           \
    void baseType<Targ>::destroy##argNames##ConstructorMaps()                  \
    {                                                                          \
        if (baseType<Targ>::argNames##ConstructorMapPtr_)                      \
        {                                                                      \
            delete baseType<Targ>::argNames##ConstructorMapPtr_;               \
            baseType<Targ>::argNames##ConstructorMapPtr_ = nullptr;            \
        }                                                                      \
    }


//- Create pointer to hash-map of functions
//  use when baseType requires the Targ template argument
#define defineTemplatedRunTimeSelectionMapPtr(baseType,argNames,Targ)          \
                                                                               \
    /* Define the constructor function map */                                  \
    baseType<Targ>::argNames##ConstructorMap*                                  \
        baseType<Targ>::argNames##ConstructorMapPtr_ = nullptr


//- Define run-time selection map for template classes
//  use when baseType requires the Targ template argument
#define defineTemplatedRunTimeSelectionMap(baseType,argNames,Targ)             \
                                                                               \
    template<>                                                                 \
    defineTemplatedRunTimeSelectionMapPtr(baseType,argNames,Targ);             \
    template<>                                                                 \
    defineTemplatedRunTimeSelectionMapConstructor(baseType,argNames,Targ);     \
    template<>                                                                 \
    defineTemplatedRunTimeSelectionMapDestructor(baseType,argNames,Targ)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

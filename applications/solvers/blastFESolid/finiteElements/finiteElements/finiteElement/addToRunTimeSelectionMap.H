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
    Macros for easy insertion into run-time selection maps

\*---------------------------------------------------------------------------*/

#ifndef addToRunTimeSelectionMap_H
#define addToRunTimeSelectionMap_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Add to hash-map of functions with typename as the key
#define addToRunTimeSelectionMap\
(baseType,thisType,argNames,lookup)                                            \
                                                                               \
    /* Add the thisType constructor function to the map */                     \
    baseType::add##argNames##ConstructorToMap<thisType>                        \
        add##thisType##argNames##ConstructorTo##baseType##Map_(lookup)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Add to hash-map of functions with typename as the key.
//  Use when baseType doesn't need a template argument (eg, is a typedef)
#define addTemplateToRunTimeSelectionMap\
(baseType,thisType,Targ,argNames,lookup)                                       \
                                                                               \
    /* Add the thisType constructor function to the map */                     \
    baseType::add##argNames##ConstructorToMap<thisType<Targ>>                  \
        add##thisType##Targ##argNames##ConstructorTo##baseType##Map_(lookup)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Add to hash-map of functions with typename as the key.
//  Use when baseType requires the Targ template argument as well
#define addTemplatedToRunTimeSelectionMap\
(baseType,thisType,Targ,argNames,lookup)                                       \
                                                                               \
    /* Add the thisType constructor function to the map */                     \
    baseType<Targ>::add##argNames##ConstructorToMap<thisType<Targ>>            \
        add##thisType##Targ##argNames##ConstructorTo##baseType##               \
        Targ##Map_(lookup)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

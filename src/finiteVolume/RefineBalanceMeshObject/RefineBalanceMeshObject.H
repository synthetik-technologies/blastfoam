/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-12-11 Jeff Heylmun:    Added storage of mixture fields and derivation
                            from UpdateableMeshObject
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

Class
    Foam::BalanceMeshObject
    Foam::RefineMeshObject
    Foam::BalanceMeshObject

Description
    Object used to lookup after refinement and balancing

SourceFiles
    RefineBalanceMeshObject.C

\*---------------------------------------------------------------------------*/

#ifndef RefineBalanceMeshObject_H
#define RefineBalanceMeshObject_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvMesh.H"
#include "regIOobject.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class RefineMeshObject Declaration
\*---------------------------------------------------------------------------*/

class BalanceMeshObject
:
    public regIOobject
{

public:

    // Declare name of the class and its debug switch
    ClassName("BalanceMeshObject");

    // Constructor
    BalanceMeshObject
    (
        const word& typeName,
        const objectRegistry& obr,
        const IOobject::readOption read = IOobject::NO_READ,
        const IOobject::writeOption write = IOobject::NO_WRITE
    )
    :
        regIOobject
        (
            IOobject
            (
                typeName,
                obr.instance(),
                obr,
                read,
                write
            )
        )
    {}

    virtual void updateObject() = 0;

    static void updateObjects(objectRegistry&);
};


class RefineMeshObject
:
    public BalanceMeshObject
{

public:

    // Declare name of the class and its debug switch
    ClassName("RefineMeshObject");

    // Constructor
    RefineMeshObject
    (
        const word& typeName,
        const objectRegistry& obr,
        const IOobject::readOption read = IOobject::NO_READ,
        const IOobject::writeOption write = IOobject::NO_WRITE
    )
    :
        BalanceMeshObject
        (
            typeName,
            obr,
            read,
            write
        )
    {}

    virtual void updateObject() = 0;

    static void updateObjects(objectRegistry&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

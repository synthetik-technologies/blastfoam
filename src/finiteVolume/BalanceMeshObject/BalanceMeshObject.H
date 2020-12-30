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

Description
    Object used to lookup after balancing and resize

SourceFiles
    BalanceMeshObject.C

\*---------------------------------------------------------------------------*/

#ifndef BalanceMeshObject_H
#define BalanceMeshObject_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvMesh.H"
#include "regIOobject.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BalanceMeshObject Declaration
\*---------------------------------------------------------------------------*/

class BalanceMeshObject
:
    public regIOobject
{

public:

    // Declare name of the class and its debug switch
    ClassName("BalanceMeshObject");

    // Constructor
    BalanceMeshObject(const word& typeName, const objectRegistry& obr)
    :
        regIOobject
        (
            IOobject
            (
                typeName,
                obr.instance(),
                obr
            )
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

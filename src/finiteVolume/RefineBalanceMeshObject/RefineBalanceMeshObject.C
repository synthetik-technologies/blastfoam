/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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

#include "RefineBalanceMeshObject.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(RefineMeshObject, 0);
    defineTypeNameAndDebug(BalanceMeshObject, 0);
}

// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

void Foam::RefineMeshObject::updateObjects(objectRegistry& obr)
{
    HashTable<RefineMeshObject*> meshObjects
    (
        obr.lookupClass<RefineMeshObject>()
    );

    if (debug)
    {
        Pout<< "RefineMeshObject::updateMesh(objectRegistry&): updating" << nl
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter
    (
        typename HashTable<RefineMeshObject*>,
        meshObjects,
        iter
    )
    {
        if (isA<RefineMeshObject>(*iter()))
        {
            if (debug)
            {
                Pout<< "    Updating " << iter()->name() << endl;
            }
            dynamic_cast<RefineMeshObject*>(iter())->updateObject();
        }
    }
}


void Foam::BalanceMeshObject::updateObjects(objectRegistry& obr)
{
    HashTable<BalanceMeshObject*> meshObjects
    (
        obr.lookupClass<BalanceMeshObject>()
    );

    if (debug)
    {
        Pout<< "BalanceMeshObject::updateMesh(objectRegistry&): updating" << nl
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter
    (
        typename HashTable<BalanceMeshObject*>,
        meshObjects,
        iter
    )
    {
        if (isA<BalanceMeshObject>(*iter()))
        {
            if (debug)
            {
                Pout<< "    Updating " << iter()->name() << endl;
            }
            dynamic_cast<BalanceMeshObject*>(iter())->updateObject();
        }
    }
}

// ************************************************************************* //

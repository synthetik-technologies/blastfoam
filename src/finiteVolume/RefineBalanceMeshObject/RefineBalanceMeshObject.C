/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020
     \\/     M anipulation  | Synthetik Applied Technology
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

#include "RefineBalanceMeshObject.H"
#include "meshObjects.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

template<class Mesh>
void Foam::blastMeshObject::preDistribute
(
    objectRegistry& obr
)
{
    HashTable<preDistributeableMeshObject<Mesh>*> meshObjects
    (
        obr.lookupClass<preDistributeableMeshObject<Mesh>>()
    );

    if (meshObjects::debug)
    {
        Pout<< "meshObject::preDistribute(objectRegistry&,"
            << "mapDistributePolyMesh&): updating " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter
    (
        typename HashTable<preDistributeableMeshObject<Mesh>*>,
        meshObjects,
        iter
    )
    {
        if (isA<preDistributeableMeshObject<Mesh>>(*iter()))
        {
            // if (meshObjects::debug)
            // {
            //     Pout<< "    Updating " << iter()->name() << endl;
            // }
            dynamic_cast<preDistributeableMeshObject<Mesh>*>
            (
                iter()
            )->preDistribute();
        }
    }
}

// ************************************************************************* //

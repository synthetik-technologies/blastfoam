/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "blastMeshTools.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "globalMeshData.H"
#include "contiguous.H"
#include "transform.H"
#include "IOobjectList.H"
#include "fvMesh.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Read and add fields to the database
template<class Type, class Mesh>
void Foam::meshTools::readInternalFields
(
    const fvMesh& mesh,
    const IOobjectList& objects
)
{
    typedef DimensionedField<Type, Mesh> dimField;
    IOobjectList fields = objects.lookupClass(dimField::typeName);
    forAllIter(IOobjectList, fields, fieldIter)
    {
        if (!mesh.foundObject<dimField>(fieldIter()->name()))
        {
            typeIOobject<dimField> fieldTargetIOobject
            (
                fieldIter()->name(),
                mesh.time().name(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            );

            if (fieldTargetIOobject.headerOk())
            {
                dimField* fPtr
                (
                    new dimField
                    (
                        fieldTargetIOobject,
                        getGeoMesh<Mesh>(mesh)
                    )
                );
                fPtr->store(fPtr);
            }
        }
    }
}


//- Read and add fields to the database
template<class Type, template<class> class Patch, class Mesh>
void Foam::meshTools::readGeoFields
(
    const fvMesh& mesh,
    const IOobjectList& objects
)
{
    typedef GeometricField<Type, Patch, Mesh> geoField;
    IOobjectList fields = objects.lookupClass(geoField::typeName);
    forAllIter(IOobjectList, fields, fieldIter)
    {
        if (!mesh.foundObject<geoField>(fieldIter()->name()))
        {
            typeIOobject<geoField> fieldTargetIOobject
            (
                fieldIter()->name(),
                mesh.time().name(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            );

            if (fieldTargetIOobject.headerOk())
            {
                geoField* fPtr
                (
                    new geoField
                    (
                        fieldTargetIOobject,
                        getGeoMesh<Mesh>(mesh)
                    )
                );
                fPtr->store(fPtr);
            }
        }
    }
}


template<class Mesh>
const typename Mesh::Mesh& Foam::meshTools::getGeoMesh(const fvMesh& mesh)
{
    return Mesh(mesh)();
}

// ************************************************************************* //

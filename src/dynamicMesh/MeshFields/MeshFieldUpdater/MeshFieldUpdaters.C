/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "MeshFieldUpdaters.H"
#include "MeshField.H"
#include "cellMapper.H"
#include "faceMapper.H"
#include "pointMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebug(cellMeshFieldUpdater, 0);

template<>
void cellMeshFieldUpdater::updateMesh(const mapPolyMesh& mpm)
{
    cellMapper mapper(mpm);
    FOR_ALL_FIELD_TYPES(updateFields, cellGeoMesh);
}

template<>
void cellMeshFieldUpdater::distribute(const mapDistributePolyMesh& mdpm)
{
    FOR_ALL_FIELD_TYPES(distributeFields, cellGeoMesh, distributeCellData);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebug(faceMeshFieldUpdater, 0);

template<>
void faceMeshFieldUpdater::updateMesh(const mapPolyMesh& mpm)
{
    faceMapper mapper(mpm);
    FOR_ALL_FIELD_TYPES(updateFields, faceGeoMesh);
}

template<>
void faceMeshFieldUpdater::distribute(const mapDistributePolyMesh& mdpm)
{
    FOR_ALL_FIELD_TYPES(distributeFields, faceGeoMesh, distributeFaceData);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebug(pointMeshFieldUpdater, 0);

template<>
void pointMeshFieldUpdater::updateMesh(const mapPolyMesh& mpm)
{
    pointMapper mapper(pointMesh::New(mesh_), mpm);
    FOR_ALL_FIELD_TYPES(updateFields, pointGeoMesh);
}

template<>
void pointMeshFieldUpdater::distribute(const mapDistributePolyMesh& mdpm)
{
    FOR_ALL_FIELD_TYPES(distributeFields, pointGeoMesh, distributePointData);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2022
     \\/     M anipulation  | Synthetik Applied Technologies
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

Description
    Iteratively set fields and refine the mesh based on a given criteria. The
    default is to check for any differences across a face, but errorEstimators
    can also be used. Selected sets can also be used to determine refinement
    zones.

    In addition to uniform values, fields can also be set using the runTime
    selectable FieldSetTypes which can be used to set non-uniform fields.
    There is no restriction on the field types that can be set
    (i.e vol/surface/point and scalar/vector/tensor/...)

    Zones and sets can also be created an modified from the region entries.

\*---------------------------------------------------------------------------*/

#include "fieldSetList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::HashTable<Foam::fieldSetList::GeoType>
Foam::fieldSetList::geoTypeNames
(
    {
        {"vol", VOL},
        {"surface", SURFACE},
        {"point", POINT}
    }
);
const Foam::HashTable
<
    Foam::word,
    Foam::fieldSetList::GeoType,
    Foam::Hash<Foam::label>
> Foam::fieldSetList::geoEnumTypes
(
    {
        {VOL, "vol"},
        {SURFACE, "surface"},
        {POINT, "point"}
    }
);

const Foam::HashTable<Foam::fieldSetList::PrimitiveType>
Foam::fieldSetList::primitiveTypeNames
(
    {
        {"Scalar", SCALAR},
        {"Vector", VECTOR},
        {"SymmTensor", SYMMTENSOR},
        {"SphericalTensor", SPHERICALTENSOR},
        {"Tensor", TENSOR}
    }
);
const Foam::HashTable
<
    Foam::word,
    Foam::fieldSetList::PrimitiveType,
    Foam::Hash<Foam::label>
> Foam::fieldSetList::primitiveEnumTypes
(
    {
        {SCALAR, "vol"},
        {VECTOR, "surface"},
        {SYMMTENSOR, "SymmTensor"},
        {SPHERICALTENSOR, "SphericalTensor"},
        {TENSOR, "Tensor"}
    }
);

Foam::fieldSetList::GeoType
Foam::fieldSetList::getGeoType(const word& type)
{
    forAllConstIter(HashTable<GeoType>, geoTypeNames, iter)
    {
        if (label(type.find(iter.key())) >= 0)
        {
            return iter();
        }
    }
    return UNKNOWN_GEO;
}


Foam::fieldSetList::PrimitiveType
Foam::fieldSetList::getPrimitiveType(const word& type)
{
    forAllConstIter(HashTable<PrimitiveType>, primitiveTypeNames, iter)
    {
        if (label(type.find(iter.key())) >= 0)
        {
            return iter();
        }
    }
    return UNKNOWN_PRIM;
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::fieldSetList::fieldSetList()
{}


Foam::fieldSetList::iNew::iNew
(
    const fvMesh& mesh,
    const dictionary& dict,
    const bool force
)
:
    mesh_(mesh),
    dict_(dict),
    write_(false),
    force_(force)
{}

Foam::fieldSetList::iNew::iNew
(
    const fvMesh& mesh,
    const dictionary& dict,
    const labelList& selectedCells,
    const labelList& selectedFaces,
    const labelList& selectedPoints,
    const bool write,
    const bool force
)
:
    mesh_(mesh),
    dict_(dict),
    selectedCells_(selectedCells),
    selectedFaces_(selectedFaces),
    selectedPoints_(selectedPoints),
    write_(write),
    force_(force)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::autoPtr<Foam::fieldSetList>
Foam::fieldSetList::iNew::operator()(Istream& is) const
{
    word fieldSetType(is);
    GeoType geo(getGeoType(fieldSetType));
    switch (geo)
    {
        case VOL:
            createTopoSet<fvPatchField, volMesh>
            (
                fieldSetType,
                selectedCells_,
                is
            );
            break;
        case SURFACE:
            createTopoSet<fvsPatchField, surfaceMesh>
            (
                fieldSetType,
                selectedFaces_,
                is
            );
            break;
        case POINT:
            createTopoSet<pointPatchField, pointMesh>
            (
                fieldSetType,
                selectedPoints_,
                is
            );
            break;
        default:
            FatalIOErrorInFunction(is)
                << "Could not determine geometry type from "
                << fieldSetType << endl
                << abort(FatalIOError);
            break;
    }

    return autoPtr<fieldSetList>(new fieldSetList());
}

// ************************************************************************* //

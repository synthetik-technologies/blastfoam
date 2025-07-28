/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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

Class
    Foam::FieldSetType

Description
    Base class used to select how a given field value is computed

SourceFiles
    fieldSetList.C

\*---------------------------------------------------------------------------*/
#ifndef fieldSetList_H
#define fieldSetList_H

#include "fvMesh.H"
#include "backupTopoSetSource.H"
#include "topoSetList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "FieldSetType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class fieldSetList Declaration
\*---------------------------------------------------------------------------*/

class fieldSetList
{

public:

    //- Geometric Type
    enum GeoType
    {
        VOL,
        SURFACE,
        POINT,
        UNKNOWN_GEO
    };
    static const HashTable<GeoType> geoTypeNames;
    static const HashTable<word, GeoType, Hash<label>> geoEnumTypes;
    static GeoType getGeoType(const word& type);


    //- Primitive Type
    enum PrimitiveType
    {
        SCALAR,
        VECTOR,
        SPHERICALTENSOR,
        SYMMTENSOR,
        TENSOR,
        UNKNOWN_PRIM
    };
    static const HashTable<PrimitiveType> primitiveTypeNames;
    static const HashTable<word, PrimitiveType, Hash<label>> primitiveEnumTypes;
    static PrimitiveType getPrimitiveType(const word& type);

public:

    //- Constructor
    fieldSetList();

    //- Clone function (required for use with PtrList)
    autoPtr<fieldSetList> clone() const
    {
        return autoPtr<fieldSetList>(new fieldSetList());
    }

    class iNew
    {
        const fvMesh& mesh_;
        const dictionary& dict_;
        const labelList selectedCells_;
        const labelList selectedFaces_;
        const labelList selectedPoints_;
        const bool write_;
        const bool force_;

        template<class Type, template<class> class Patch, class Mesh>
        void createTopoSetType
        (
            const word& fieldSetTypeDesc,
            const labelList& elms,
            Istream& is
        ) const
        {
            const word fieldName(is);
            autoPtr<FieldSetType<Type, Patch, Mesh>> fieldSet
            (
                FieldSetType<Type, Patch, Mesh>::New
                (
                    fieldSetTypeDesc,
                    fieldName,
                    mesh_,
                    dict_,
                    elms,
                    is,
                    write_
                )
            );
            if (fieldSet->good())
            {
                Info<< "    Setting " << fieldName << endl;
            }
            else
            {
                WarningInFunction
                    << "Field " << fieldName << " not found" << endl;
            }
        }

        template<template<class> class Patch, class Mesh>
        void createTopoSet
        (
            const word& fieldSetType,
            const labelList& elms,
            Istream& is
        ) const
        {
            PrimitiveType prim(getPrimitiveType(fieldSetType));
            switch (prim)
            {
                case SCALAR:
                    createTopoSetType<scalar, Patch, Mesh>
                    (
                        fieldSetType,
                        elms,
                        is
                    );
                    break;
                case VECTOR:
                    createTopoSetType<vector, Patch, Mesh>
                    (
                        fieldSetType,
                        elms,
                        is
                    );
                    break;
                case SYMMTENSOR:
                    createTopoSetType<symmTensor, Patch, Mesh>
                    (
                        fieldSetType,
                        elms,
                        is
                    );
                    break;
                case SPHERICALTENSOR:
                    createTopoSetType<sphericalTensor, Patch, Mesh>
                    (
                        fieldSetType,
                        elms,
                        is
                    );
                    break;
                case TENSOR:
                    createTopoSetType<tensor, Patch, Mesh>
                    (
                        fieldSetType,
                        elms,
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
        }


    public:

        iNew
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const bool force = true
        );

        iNew
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const labelList& selectedCells,
            const labelList& selectedFaces,
            const labelList& selectedPoints,
            const bool write,
            const bool force = false
        );

        autoPtr<fieldSetList> operator()(Istream& is) const;
    };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

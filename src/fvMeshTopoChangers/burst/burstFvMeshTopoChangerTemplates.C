/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
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

\*---------------------------------------------------------------------------*/

#include "burstFvMeshTopoChanger.H"
#include "conformedFvPatchField.H"
#include "conformedFvsPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
void Foam::fvMeshTopoChangers::burst::storeVolBoundaries
(
    const labelList& boundaryMap,
    HashPtrTable<typename VolField<Type>::Boundary>& bfields
) const
{
    const fvBoundaryMesh& bMesh = mesh().boundary();
    UPtrList<VolField<Type>> fields(mesh().curFields<VolField<Type>>());
    forAll(fields, i)
    {
        VolField<Type>& field = fields[i];
        bfields.insert
        (
            field.name(),
            new typename VolField<Type>::Boundary(bMesh)
        );
        typename VolField<Type>::Boundary& bfield0 = *bfields[field.name()];
        typename VolField<Type>::Boundary& bfield =
            const_cast<typename VolField<Type>::Boundary&>(field.boundaryField());

        // Unconform mesh, i.e. make sure the actual patch types are being used
        // since the old patch is not conformed but the new patch is
        conformedFvPatchField<Type>::unconform(bfield);

        forAll(boundaryMap, patchi)
        {
            const label newPatchi = boundaryMap[patchi];
            if (newPatchi >= 0)
            {
                bfield0.set
                (
                    newPatchi,
                    field.boundaryField()[patchi].clone(field)
                );
            }
        }
    }
}


template<class Type>
void Foam::fvMeshTopoChangers::burst::storeSurfaceBoundaries
(
    const labelList& boundaryMap,
    HashPtrTable<typename SurfaceField<Type>::Boundary>& bfields
) const
{
    const fvBoundaryMesh& bMesh = mesh().boundary();
    UPtrList<SurfaceField<Type>> fields(mesh().curFields<SurfaceField<Type>>());
    forAll(fields, i)
    {
        SurfaceField<Type>& field = fields[i];
        bfields.insert
        (
            field.name(),
            new typename SurfaceField<Type>::Boundary(bMesh)
        );
        typename SurfaceField<Type>::Boundary& bfield0 = *bfields[field.name()];
        typename SurfaceField<Type>::Boundary& bfield =
            const_cast<typename SurfaceField<Type>::Boundary&>(field.boundaryField());

        // Unconform mesh, i.e. make sure the actual patch types are being used
        // since the old patch is not conformed but the new patch is
        conformedFvsPatchField<Type>::unconform(bfield);

        forAll(boundaryMap, patchi)
        {
            const label newPatchi = boundaryMap[patchi];
            if (newPatchi >= 0)
            {
                bfield0.set
                (
                    newPatchi,
                    field.boundaryField()[patchi].clone(field)
                );
            }
        }
    }
}


template<class Type>
void Foam::fvMeshTopoChangers::burst::mapVolBoundaries
(
    const PtrList<fieldMapper>& mappers,
    const HashPtrTable<typename VolField<Type>::Boundary>& bfields
) const
{
    UPtrList<VolField<Type>> fields(mesh().curFields<VolField<Type>>());
    forAll(fields, i)
    {
        VolField<Type>& field = fields[i];
        typename VolField<Type>::Boundary& bfield =
            const_cast<typename VolField<Type>::Boundary&>(field.boundaryField());
        const typename VolField<Type>::Boundary& bfield0 = *bfields[field.name()];

        forAll(bfield0, patchi)
        {
            if (bfield0.set(patchi))
            {
                fvPatchField<Type>& fvp = const_cast<fvPatchField<Type>&>
                (
                    field.boundaryField()[patchi]
                );
                fvp.map(bfield0[patchi], mappers[patchi]);
            }
        }

        // Mapping has been completed so re-conform the patch fields
        conformedFvPatchField<Type>::conform(bfield);
    }
}


template<class Type>
void Foam::fvMeshTopoChangers::burst::mapSurfaceBoundaries
(
    const PtrList<fieldMapper>& mappers,
    const HashPtrTable<typename SurfaceField<Type>::Boundary>& bfields
) const
{
    UPtrList<SurfaceField<Type>> fields(mesh().curFields<SurfaceField<Type>>());
    forAll(fields, i)
    {
        SurfaceField<Type>& field = fields[i];
        typename SurfaceField<Type>::Boundary& bfield =
            const_cast<typename SurfaceField<Type>::Boundary&>(field.boundaryField());
        const typename SurfaceField<Type>::Boundary& bfield0 = *bfields[field.name()];

        forAll(bfield0, patchi)
        {
            if (bfield0.set(patchi))
            {
                fvsPatchField<Type>& fvsp = const_cast<fvsPatchField<Type>&>
                (
                    field.boundaryField()[patchi]
                );
                fvsp.map(bfield0[patchi], mappers[patchi]);
            }
        }

        // Mapping has been completed so re-conform the patch fields
        conformedFvsPatchField<Type>::conform(bfield);
    }
}


template<class Type>
void Foam::fvMeshTopoChangers::burst::setUnmappedValues
(
    const PackedBoolList& mappedFace
) const
{
    UPtrList<VolField<Type>> fields(mesh().curFields<VolField<Type>>());

    forAll(fields, i)
    {
        VolField<Type>& field = fields[i];

        forAll(field.boundaryField(), patchi)
        {
            fvPatchField<Type>& fvp = const_cast<fvPatchField<Type>&>
            (
                field.boundaryField()[patchi]
            );
            const label start = fvp.patch().start();
            const labelList& faceCells = fvp.patch().faceCells();
            forAll(fvp, fi)
            {
                if (!mappedFace[start+fi])
                {
                    fvp[fi] = field[faceCells[fi]];
                }
            }
        }
    }
}

// ************************************************************************* //

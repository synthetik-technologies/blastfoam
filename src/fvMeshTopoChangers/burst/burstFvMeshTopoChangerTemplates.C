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

#include "fixedGradientFvPatchField.H"
#include "mixedFvPatchField.H"
#include "directionMixedFvPatchField.H"

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

        // Unconform patch field, i.e. make sure the actual patch types are
        // being used since the old patch is not conformed but the new patch is
        conformedFvPatchField<Type>::unconform(bfield);

        // Store
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

        // Unconform patch field, i.e. make sure the actual patch types are
        // being used since the old patch is not conformed but the new patch is
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
    const List<labelList>& addressing,
    const HashPtrTable<typename VolField<Type>::Boundary>& bfields
) const
{
    auto evaluate = [](const typename VolField<Type>::Patch& pf)
    {
        return
            (
                isA<nonConformalFvPatch>(pf.patch())
             && pf.type() == pf.patch().patch().type()
             && polyPatch::constraintType(pf.patch().patch().type())
            )
         || isA<nonConformalErrorFvPatch>(pf.patch());
    };

    UPtrList<VolField<Type>> fields(mesh().curFields<VolField<Type>>());
    forAll(fields, i)
    {
        VolField<Type>& field = fields[i];
        typename VolField<Type>::Boundary& bfield =
            field.boundaryFieldRefNoStoreOldTimes();
        const typename VolField<Type>::Boundary& bfield0 =
            *bfields[field.name()];

        forAll(bfield0, patchi)
        {
            const labelList& addr = addressing[patchi];
            if (addr.size())
            {
                forwardFieldMapper m(addr);

                typename VolField<Type>::Patch& pf = bfield[patchi];
                const typename VolField<Type>::Patch& pf0 = bfield0[patchi];

                // Check for comparable types between new and old patch fields
                // If the types are not the same then check basic types
                if
                (
                    (pf.type() == pf0.type())
                 || (
                        isA<fixedGradientFvPatchField<Type>>(pf)
                     && isA<fixedGradientFvPatchField<Type>>(pf0)
                    )
                 || (
                        isA<mixedFvPatchField<Type>>(pf)
                     && isA<mixedFvPatchField<Type>>(pf0)
                    )
                 || (
                        isA<directionMixedFvPatchField<Type>>(pf)
                     && isA<directionMixedFvPatchField<Type>>(pf0)
                    )
                )
                {
                    // Only map if same type
                    pf.map(pf0, m);
                }
                else
                {
                    // Map the actual values
                    m(pf, pf0);

                    // Check type to try and determine additional fields to
                    // map since actual values will be unmapped
                    if (isA<fixedGradientFvPatchField<Type>>(pf))
                    {
                        fixedGradientFvPatchField<Type>& fgpf =
                            dynamicCast<fixedGradientFvPatchField<Type>>(pf);
                        Field<Type>& g = fgpf.gradient();
                        forAll(addr, fi)
                        {
                            if (addr[fi] >= 0)
                            {
                                g[fi] = Zero;
                            }
                        }
                    }
                    else if (isA<mixedFvPatchField<Type>>(pf))
                    {
                        mixedFvPatchField<Type>& mpf =
                            dynamicCast<mixedFvPatchField<Type>>(pf);
                        Field<Type>& rv = mpf.refValue();
                        Field<Type>& rg = mpf.refGrad();
                        Field<scalar>& vf = mpf.valueFraction();
                        forAll(addr, fi)
                        {
                            if (addr[fi] >= 0)
                            {
                                // Set ref value as the current value
                                rv[fi] = mpf[fi];

                                // Zero gradient
                                rg[fi] = Zero;
                                vf[fi] = Zero;
                            }
                        }
                    }
                    else if (isA<directionMixedFvPatchField<Type>>(pf))
                    {
                        directionMixedFvPatchField<Type>& dmpf =
                            dynamicCast<directionMixedFvPatchField<Type>>(pf);
                        Field<Type>& rv = dmpf.refValue();
                        Field<Type>& rg = dmpf.refGrad();
                        Field<symmTensor>& vf = dmpf.valueFraction();
                        forAll(addr, fi)
                        {
                            if (addr[fi] >= 0)
                            {
                                // Set ref value as the current value
                                rv[fi] = dmpf[fi];

                                // Zero gradient
                                rg[fi] = Zero;
                                vf[fi] = Zero;
                            }
                        }
                    }
                    // Not a basic boundary so not sure how to fix
                    // Un-initialize fields. Hopefully this is okay
                }
            }
        }

        // Mapping has been completed so re-conform the patch fields
        conformedFvPatchField<Type>::conform(bfield);


        // Synchronise boundaries
        const label nReq = Pstream::nRequests();

        forAll(bfield0, patchi)
        {
            typename VolField<Type>::Patch& pf = bfield[patchi];
            if (bfield0.set(patchi) && evaluate(pf))
            {
                pf.initEvaluate(Pstream::defaultCommsType);
            }
        }

        if
        (
            Pstream::parRun()
         && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            Pstream::waitRequests(nReq);
        }

        forAll(bfield0, patchi)
        {
            typename VolField<Type>::Patch& pf = bfield[patchi];
            if (bfield0.set(patchi) && evaluate(pf))
            {
                pf.evaluate(Pstream::defaultCommsType);
            }
        }
    }
}


template<class Type>
void Foam::fvMeshTopoChangers::burst::mapSurfaceBoundaries
(
    const List<labelList>& addressing,
    const HashPtrTable<typename SurfaceField<Type>::Boundary>& bfields
) const
{
    UPtrList<SurfaceField<Type>> fields(mesh().curFields<SurfaceField<Type>>());
    forAll(fields, i)
    {
        SurfaceField<Type>& field = fields[i];
        typename SurfaceField<Type>::Boundary& bfield =
            field.boundaryFieldRefNoStoreOldTimes();
        const typename SurfaceField<Type>::Boundary& bfield0 =
            *bfields[field.name()];

        forAll(bfield0, patchi)
        {
            const labelList& addr = addressing[patchi];
            if (addr.size())
            {
                typename SurfaceField<Type>::Patch& pf = bfield[patchi];
                pf.map(bfield0[patchi], forwardFieldMapper(addr));
            }
        }

        // Mapping has been completed so re-conform the patch fields
        conformedFvsPatchField<Type>::conform(bfield);

        // bfield = fvMeshStitcherTools::synchronisedBoundaryField
        // (
        //     bfield,
        //     false,
        //     0.5,
        //     0.5
        // );
    }
}

// ************************************************************************* //

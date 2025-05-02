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

#include "localLeastSquaresGrad.H"
#include "localLeastSquaresVectors.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
localLeastSquaresGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vsf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tlsGrad
    (
        GeometricField<GradType, fvPatchField, volMesh>::New
        (
            name,
            mesh,
            dimensioned<GradType>
            (
                "zero",
                vsf.dimensions()/dimLength,
                pTraits<GradType>::zero
            )
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& lsGrad =
        tlsGrad.ref();

    // Get reference to least square vectors
    const localLeastSquaresVectors& lsv =
        localLeastSquaresVectors::New
        (
            mesh,
            local_
        );

    const tensorField& Ainv = lsv.Ainv();

    const volVectorField& C = mesh.C();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        GradType deltaVsf((vsf[nei] - vsf[own])*(C[nei] - C[own]));

        lsGrad[own] += Ainv[own] & deltaVsf;
        lsGrad[nei] += Ainv[nei] & deltaVsf;
    }

    // Boundary faces
    typename GeometricField<GradType, fvPatchField, volMesh>::Boundary& blsGrad =
        lsGrad.boundaryFieldRef();
    forAll(blsGrad, patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];

        if (vsf.boundaryField()[patchi].coupled())
        {
            const vectorField pdelta(patch.delta());
            const labelList& faceCells = patch.faceCells();
            Field<Type> neiVsf
            (
                vsf.boundaryField()[patchi].patchNeighbourField()
            );

            // blsGrad[patchi] =
            //     patch.deltaCoeffs()
            //    *(neiVsf - pU[patchi].patchInternalField())*patch.nf();

            forAll(patch, facei)
            {
                const label& celli = faceCells[facei];

                lsGrad[celli] +=
                    Ainv[celli] & ((neiVsf[facei] - vsf[celli])*pdelta[facei]);
            }
        }
        // else
        // {
        //     const fvPatchField<Type>& patchVsf = vsf.boundaryField()[patchi];
        //
        //     forAll(patchVsf, patchFaceI)
        //     {
        //         lsGrad[faceCells[patchFaceI]] +=
        //              patchOwnLs[patchFaceI]
        //             *(patchVsf[patchFaceI] - vsf[faceCells[patchFaceI]]);
        //     }
        // }
    }

    // lsGrad.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vsf, lsGrad);

    return tlsGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

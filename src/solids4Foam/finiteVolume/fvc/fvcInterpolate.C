/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "fvcGradf.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "ggiFvPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField<Type, pointPatchField, pointMesh>& pf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tResult
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "interpolate(" + pf.name() + ")",
                pf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>
            (
                "0",
                vf.dimensions(),
                pTraits<Type>::zero
            )
        )
    );

    Field<Type>& resultI = tResult().internalField();


    const vectorField& points = mesh.points();

    const faceList& faces = mesh.faces();

    const Field<Type>& pfI = pf.internalField();

//     const unallocLabelList& owner = mesh.owner();
//     const unallocLabelList& neighbour = mesh.neighbour();

    forAll(resultI, faceI)
    {
        const face& curFace = faces[faceI];

        // If the face is a triangle, do a direct calculation
        if (curFace.size() == 3)
        {
            resultI[faceI] = curFace.average(points, pfI);
        }
        else
        {
            label nPoints = curFace.size();

            point centrePoint = point::zero;
            Type cf = pTraits<Type>::zero;

            for (register label pI=0; pI<nPoints; pI++)
            {
                centrePoint += points[curFace[pI]];
                cf += pfI[curFace[pI]];
            }

            centrePoint /= nPoints;
            cf /= nPoints;

            resultI[faceI] = cf;
        }
    }

    forAll(mesh.boundary(), patchI)
    {
        tResult().boundaryField()[patchI] =
            vf.boundaryField()[patchI];

//         forAll(mesh.boundary()[patchI], faceI)
//         {
//             label globalFaceID =
//                 mesh.boundaryMesh()[patchI].start() + faceI;

//             const face& curFace = faces[globalFaceID];

//             // If the face is a triangle, do a direct calculation
//             if (curFace.size() == 3)
//             {
//                 tResult().boundaryField()[patchI][faceI] =
//                     curFace.average(points, pfI);
//             }
//             else
//             {
//                 label nPoints = curFace.size();

//                 point centrePoint = point::zero;
//                 Type cf = pTraits<Type>::zero;

//                 for (register label pI=0; pI<nPoints; pI++)
//                 {
//                     centrePoint += points[curFace[pI]];
//                     cf += pfI[curFace[pI]];
//                 }

//                 centrePoint /= nPoints;
//                 cf /= nPoints;

//                 tResult().boundaryField()[patchI][faceI] = cf;
//             }
//         }
    }

    return tResult;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "globalPolyPatch.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::faceToPoint
(
    const Field<Type>& fFld
) const
{
    const primitivePatch& patch = physicalPatch();
    if (fFld.size() != patch.size())
    {
        FatalErrorInFunction
            << "Patch field does not correspond to patch faces.  Patch size: "
            << patch.size() << " field size: " << fFld.size()
            << abort(FatalError);
    }

    tmp<Field<Type>> tpFld(new Field<Type>(patch.nPoints()));
    Field<Type>& pFld = tpFld.ref();

    // point to face addressing
    const labelListList& pointFaces = patch.pointFaces();

    // Compute face to point weights (inverse distance)
    const List<scalarField>& weights = faceToPointWeights();
    forAll(pointFaces, pointi)
    {
        const labelList& pfs = pointFaces[pointi];
        const scalarField& ws = weights[pointi];
        Type sumWF = Zero;
        forAll(pfs, pfi)
        {
            sumWF += ws[pfi]*fFld[pfs[pfi]];
        }
        pFld[pointi] = sumWF;
    }

    syncTools::syncPointList
    (
        mesh_,
        polyPatch_.meshPoints(),
        pFld,
        plusEqOp<Type>(),
        pTraits<Type>::zero
    );

    pFld /= faceToPointSumWeights();
    return tpFld;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::faceToPoint
(
    const tmp<Field<Type>>& tfFld
) const
{
    return faceToPoint(tfFld());
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::pointToFace
(
    const Field<Type>& pFld
) const
{
    if (pFld.size() != polyPatch_.nPoints())
    {
        FatalErrorInFunction
            << "Patch field does not correspond to patch points. "
            << "Patch size: "
            << physicalPatch().nPoints() << " field size: " << pFld.size()
            << abort(FatalError);
    }

    return pointToFaceInterpolator().pointToFaceInterpolate(pFld);
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::pointToFace
(
    const tmp<Field<Type>>& pField
) const
{
    return pointToFace(pField());
}


// ************************************************************************* //

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::patchPointToGlobal
(
    const Field<Type>& pField
) const
{
    if (pField.size() != patch().nPoints())
    {
        FatalErrorInFunction
            << "Patch field does not correspond to patch points.  Patch size: "
            << patch().nPoints() << " field size: " << pField.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tgField
    (
        new Field<Type>(globalPatch().nPoints(), pTraits<Type>::zero)
    );
    Field<Type>& gField = tgField.ref();


    if (Pstream::parRun())
    {
        // PC, 16/12/17
        // We have removed duplicate points so multiple local processor points
        // may map to the same global point, which we will account for using
        // the nPoints field
        scalarField nPoints(gField.size(), 0.0);

        const labelList& addr = pointToGlobalAddr();

        forAll(addr, i)
        {
            const label globalPointID = addr[i];
            gField[globalPointID] = pField[i];
            nPoints[globalPointID] += 1.0;
        }

        // Global comm
        reduce(gField, sumOp<List<Type> >());
        reduce(nPoints, sumOp<List<scalar> >());
        gField /= nPoints;
    }
    else
    {
        gField = pField;
    }
    return tgField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::patchPointToGlobal
(
    const tmp<Field<Type>>& pField
) const
{
    return patchPointToGlobal(pField());
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::globalPointToPatch
(
    const Field<Type>& gField
) const
{
    if (gField.size() != globalPatch().nPoints())
    {
        FatalErrorInFunction
            << "Patch field does not correspond to global patch points.  "
            << "Global patch size: " << globalPatch().nPoints()
            << " field size: " << gField.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tpField
    (
        new Field<Type>(patch().nPoints(), pTraits<Type>::zero)
    );
    Field<Type>& pField = tpField.ref();

    if (Pstream::parRun())
    {
        const labelList& addr = pointToGlobalAddr();

        forAll (addr, i)
        {
            pField[i] = gField[addr[i]];
        }
    }
    else
    {
        pField = gField;
    }

    return tpField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::globalPointToPatch
(
    const tmp<Field<Type>>& gField
) const
{
    return globalPointToPatch(gField());
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::patchFaceToGlobal
(
    const Field<Type>& pField
) const
{
    if (pField.size() != patch().size())
    {
        FatalErrorInFunction
            << "Patch field does not correspond to patch faces.  Patch size: "
            << patch().size() << " field size: " << pField.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tgField
    (
        new Field<Type>(globalPatch().size(), pTraits<Type>::zero)
    );
    Field<Type>& gField = tgField.ref();

    if (Pstream::parRun())
    {
        const labelList& addr = faceToGlobalAddr();

        forAll (addr, i)
        {
            gField[addr[i]] = pField[i];
        }

        // Global comm
        reduce(gField, sumOp<List<Type> >());
    }
    else
    {
        gField = pField;
    }

    return tgField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::patchFaceToGlobal
(
    const tmp<Field<Type>>& pField
) const
{
    return patchFaceToGlobal(pField());
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::globalFaceToPatch
(
    const Field<Type>& gField
) const
{
    if (gField.size() != globalPatch().size())
    {
        FatalErrorInFunction
            << "Patch field does not correspond to global patch faces.  "
            << "Global patch size: " << globalPatch().size()
            << " field size: " << gField.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tpField
    (
        new Field<Type>(patch().size(), pTraits<Type>::zero)
    );
    Field<Type>& pField = tpField.ref();

    if (Pstream::parRun())
    {
        const labelList& addr = faceToGlobalAddr();

        forAll(addr, i)
        {
            pField[i] = gField[addr[i]];
        }
    }
    else
    {
        pField = gField;
    }

    return tpField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::globalFaceToPatch
(
    const tmp<Field<Type>>& gField
) const
{
    return globalFaceToPatch(gField());
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::faceToPoint
(
    const Field<Type>& fField
) const
{
    if (fField.size() != patch().size())
    {
        FatalErrorInFunction
            << "Patch field does not correspond to patch faces.  Patch size: "
            << patch().size() << " field size: " << fField.size()
            << abort(FatalError);
    }

    if (Pstream::parRun())
    {
        return
            globalPointToPatch
            (
                interpPtr_->faceToPointInterpolate
                (
                    patchFaceToGlobal(fField)
                )
            );
    }
    return interpPtr_->faceToPointInterpolate(fField);
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::faceToPoint
(
    const tmp<Field<Type>>& fField
) const
{
    return faceToPoint(fField());
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::pointToFace
(
    const Field<Type>& pField
) const
{
    if (pField.size() != patch().nPoints())
    {
        FatalErrorInFunction
            << "Patch field does not correspond to patch points.  Patch size: "
            << patch().nPoints() << " field size: " << pField.size()
            << abort(FatalError);
    }

    return localInterpolator().pointToFaceInterpolate(pField());
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

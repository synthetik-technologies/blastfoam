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

#include "newLeastSquaresVolPointInterpolation.H"
#include "volFields.H"
#include "pointFields.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "cyclicFvPatch.H"
#include "transform.H"
#ifdef FOAMEXTEND
    #include "cyclicGgiPolyPatch.H"
    #include "cyclicGgiFvPatchFields.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
void newLeastSquaresVolPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::interpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    const Field<Type>& vfI = vf.primitiveField();
    Field<Type>& pfI = pf.primitiveFieldRef();

    const vectorField& points = mesh().points();

    const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();

    label nCoeffs = 3;
    const labelListList& ptCells = mesh().pointCells();
    const labelListList& ptBndFaces = pointBndFaces();
    const labelListList& ptCyclicFaces = pointCyclicFaces();
    const labelListList& ptCyclicBndFaces = pointCyclicBndFaces();
    const labelListList& ptProcFaces = pointProcFaces();

    Map<Field<Type> > gPtNgbProcBndFaceFieldData;
    globalPointNgbProcBndFaceFieldData(vf, gPtNgbProcBndFaceFieldData);

    Map<Field<Type> > gPtNgbProcCellFieldData;
    globalPointNgbProcCellFieldData(vf, gPtNgbProcCellFieldData);

    const Map<List<labelPair> >& ptProcCells = pointProcCells();

    const List<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();
//     const Map<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();

    const FieldField<Field, scalar>& w = weights();
    const vectorField& o = origins();
    const scalarField& L = refL();

    FieldField<Field, Type> procCellVfI = procCellsFieldData(vfI);
    FieldField<Field, Type> procBndFaceVf = procBndFacesFieldData(vf);

    forAll(pfI, pointI)
    {
        const scalarRectangularMatrix& curInvMatrix = invMatrices[pointI];

        const scalarField& W = w[pointI];

        const labelList& interpCells = ptCells[pointI];
        const labelList& interpBndFaces = ptBndFaces[pointI];
        const labelList& interpCyclicFaces = ptCyclicFaces[pointI];
        const labelList& interpCyclicBndFaces = ptCyclicBndFaces[pointI];
        const labelList& interpProcFaces = ptProcFaces[pointI];

        Field<Type> glInterpNgbProcBndFaceData(0);
        if (gPtNgbProcBndFaceFieldData.found(pointI))
        {
            glInterpNgbProcBndFaceData = gPtNgbProcBndFaceFieldData[pointI];
        }

        Field<Type> glInterpNgbProcCellData(0);
        if (gPtNgbProcCellFieldData.found(pointI))
        {
            glInterpNgbProcCellData = gPtNgbProcCellFieldData[pointI];
        }

        Field<Type> interpNgbProcCellData(0);
        if (ptProcCells.found(pointI))
        {
            const List<labelPair>& pc = ptProcCells[pointI];

            interpNgbProcCellData.setSize(pc.size());

            forAll(pc, cI)
            {
                interpNgbProcCellData[cI] =
                    procCellVfI
                    [
                        pc[cI].first()
                    ]
                    [
                        pc[cI].second()
                    ];
            }
        }

        Field<Type> interpNgbProcBndFaceData(0);
        if (ptProcBndFaces[pointI].size())
//         if (ptProcBndFaces.found(pointI))
        {
            const List<labelPair>& pf = ptProcBndFaces[pointI];

            interpNgbProcBndFaceData.setSize(pf.size());

            forAll(pf, fI)
            {
                interpNgbProcBndFaceData[fI] =
                    procBndFaceVf
                    [
                        pf[fI].first()
                    ]
                    [
                        pf[fI].second()
                    ];
            }
        }

        Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
        Field<Type> source
        (
            interpCells.size()
          + interpBndFaces.size()
          + interpCyclicFaces.size()
          + interpCyclicBndFaces.size()
          + interpProcFaces.size()
          + glInterpNgbProcBndFaceData.size()
          + glInterpNgbProcCellData.size()
          + interpNgbProcCellData.size()
          + interpNgbProcBndFaceData.size(),
            pTraits<Type>::zero
        );

        label pointID = 0;

        Type avg = pTraits<Type>::zero;

        for (label i=0; i<interpCells.size(); i++)
        {
            source[pointID] = vfI[interpCells[i]];
            avg += sqr(W[pointID])*vfI[interpCells[i]];
            pointID++;
        }

        for (label i=0; i<interpBndFaces.size(); i++)
        {
            label faceID = interpBndFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID] = vf.boundaryField()[patchID][localFaceID];
            avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
            pointID++;
        }

        for (label i=0; i<interpCyclicFaces.size(); i++)
        {
            label faceID = interpCyclicFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            const unallocLabelList& faceCells =
                mesh().boundary()[patchID].faceCells();

            const cyclicFvPatch& cycPatch =
                refCast<const cyclicFvPatch>(mesh().boundary()[patchID]);

            label sizeby2 = faceCells.size()/2;

            if
            (
               !(
                    cycPatch.transform().transformsPosition()
                 || pTraits<Type>::rank == 0)
                )
            {
                if (localFaceID < sizeby2)
                {
                    source[pointID] =
                        transform
                        (
                            cycPatch,
                            cycPatch.transform().T()[0],
                            vfI[faceCells[localFaceID + sizeby2]]
                        );
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
                else
                {
                    source[pointID] =
                        cycPatch.invTransformPoints
                        (
                            vfI[faceCells[localFaceID - sizeby2]],
                            vfI[faceCells[localFaceID - sizeby2]]
                        );
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
            }
            else
            {
                if (localFaceID < sizeby2)
                {
                    source[pointID] = vfI[faceCells[localFaceID + sizeby2]];
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
                else
                {
                    source[pointID] = vfI[faceCells[localFaceID - sizeby2]];
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
            }
        }

        for (label i=0; i<interpCyclicBndFaces.size(); i++)
        {
            label cycFaceID = interpCyclicFaces[0];
            label cycPatchID = mesh().boundaryMesh().whichPatch(cycFaceID);
            label cycStart = mesh().boundaryMesh()[cycPatchID].start();
            label cycLocalFaceID = cycFaceID - cycStart;

            const cyclicPolyPatch& cycPolyPatch =
                refCast<const cyclicPolyPatch>
                (
                    mesh().boundaryMesh()[cycPatchID]
                );

            label faceID = interpCyclicBndFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);
            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            label sizeby2 = cycPolyPatch.size()/2;

            if (!cycPolyPatch.parallel())
            {
                if (cycLocalFaceID < sizeby2)
                {
                    source[pointID] =
                        transform
                        (
                            cycPolyPatch.forwardT()[0],
                            vf.boundaryField()[patchID][localFaceID]
                        );
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
                else
                {
                    source[pointID] =
                        transform
                        (
                            cycPolyPatch.reverseT()[0],
                            vf.boundaryField()[patchID][localFaceID]
                        );
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
            }
            else
            {
                if (cycLocalFaceID < sizeby2)
                {
                    source[pointID] =
                        vf.boundaryField()[patchID][localFaceID];
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
                else
                {
                    source[pointID] =
                        vf.boundaryField()[patchID][localFaceID];
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
            }
        }

        for (label i=0; i<interpProcFaces.size(); i++)
        {
            label faceID = interpProcFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID] = vf.boundaryField()[patchID][localFaceID];
            avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
            pointID++;
        }

        for (label i=0; i<glInterpNgbProcBndFaceData.size(); i++)
        {
            source[pointID] = glInterpNgbProcBndFaceData[i];
            avg += sqr(W[pointID])*glInterpNgbProcBndFaceData[i];
            pointID++;
        }

        for (label i=0; i<glInterpNgbProcCellData.size(); i++)
        {
            source[pointID] = glInterpNgbProcCellData[i];
            avg += sqr(W[pointID])*glInterpNgbProcCellData[i];
            pointID++;
        }

        for (label i=0; i<interpNgbProcCellData.size(); i++)
        {
            source[pointID] = interpNgbProcCellData[i];
            avg += sqr(W[pointID])*interpNgbProcCellData[i];
            pointID++;
        }

        for (label i=0; i<interpNgbProcBndFaceData.size(); i++)
        {
            source[pointID] = interpNgbProcBndFaceData[i];
            avg += sqr(W[pointID])*interpNgbProcBndFaceData[i];
            pointID++;
        }

        if (mag(mirrorPlaneTransformation()[pointI].first()) > SMALL)
//         if (mirrorPlaneTransformation().found(pointI))
        {
            const tensor& T = mirrorPlaneTransformation()[pointI].second();

            label oldSize = source.size();

            source.setSize(2*oldSize);

            for (label i=oldSize; i<source.size(); i++)
            {
                source[i] = transform(T, source[i-oldSize]);
            }

            avg += transform(T, avg);
        }

//         avg /= source.size() + SMALL;
        avg /= sum(sqr(W));

        source -= avg;

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

        vector dr = (points[pointI] - o[pointI])/L[pointI];

        pfI[pointI] =
            avg
          + coeffs[0]*dr.x()
          + coeffs[1]*dr.y()
          + coeffs[2]*dr.z();
    }

    pf.correctBoundaryConditions();

    // Correct axis point values
    forAll(mesh().boundaryMesh(), patchI)
    {
        if
        (
            mesh().boundaryMesh()[patchI].type()
         == wedgePolyPatch::typeName
        )
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            const wedgePolyPatch& wedge =
                refCast<const wedgePolyPatch>
                (
                    mesh().boundaryMesh()[patchI]
                );

            vector n =
                transform
                (
                    wedge.faceT(),
                    wedge.centreNormal()
                );
            n /= mag(n);

            forAll(meshPoints, pointI)
            {
                label curMeshPoint = meshPoints[pointI];

                if (pointAxisEdges().found(curMeshPoint))
                {
                    // Keep only component along axis
                    pfI[curMeshPoint] =
                        transform
                        (
                            sqr(wedge.axis()),
                            pfI[curMeshPoint]
                        );
                }
                else
                {
                    pfI[curMeshPoint] =
                        transform
                        (
                            I-sqr(n),
                            pfI[curMeshPoint]
                        );
                }
            }
        }
    }
}


template<class Type>
tmp<Field<Type> > newLeastSquaresVolPointInterpolation::interpolate
(
    const polyPatch& patch,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::interpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    const Field<Type>& vfI = vf.primitiveField();

    tmp<Field<Type> > tppf
    (
        new Field<Type>
        (
            patch.nPoints(),
            pTraits<Type>::zero
        )
    );
    Field<Type>& ppf = tppf.ref();

    const labelList& meshPoints = patch.meshPoints();

    const vectorField& points = mesh().points();

    const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();

    label nCoeffs = 3;
    const labelListList& ptCells = mesh().pointCells();
    const labelListList& ptBndFaces = pointBndFaces();
    const labelListList& ptCyclicFaces = pointCyclicFaces();
    const labelListList& ptProcFaces = pointProcFaces();

    Map<Field<Type> > gPtNgbProcBndFaceFieldData;
    globalPointNgbProcBndFaceFieldData(vf, gPtNgbProcBndFaceFieldData);

    Map<Field<Type> > gPtNgbProcCellFieldData;
    globalPointNgbProcCellFieldData(vf, gPtNgbProcCellFieldData);

    const Map<List<labelPair> >& ptProcCells = pointProcCells();

    const List<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();
//     const Map<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();

    const FieldField<Field, scalar>& w = weights();
    const vectorField& o = origins();
    const scalarField& L = refL();

    FieldField<Field, Type> procCellVfI = procCellsFieldData(vfI);
    FieldField<Field, Type> procBndFaceVf = procBndFacesFieldData(vf);

    forAll(ppf, pI)
    {
        label pointI = meshPoints[pI];

        const scalarRectangularMatrix& curInvMatrix = invMatrices[pointI];

        const scalarField& W = w[pointI];

        const labelList& interpCells = ptCells[pointI];
        const labelList& interpBndFaces = ptBndFaces[pointI];
        const labelList& interpCyclicFaces = ptCyclicFaces[pointI];
        const labelList& interpProcFaces = ptProcFaces[pointI];

        Field<Type> glInterpNgbProcBndFaceData(0);
        if (gPtNgbProcBndFaceFieldData.found(pointI))
        {
            glInterpNgbProcBndFaceData = gPtNgbProcBndFaceFieldData[pointI];
        }

        Field<Type> glInterpNgbProcCellData(0);
        if (gPtNgbProcCellFieldData.found(pointI))
        {
            glInterpNgbProcCellData = gPtNgbProcCellFieldData[pointI];
        }

        Field<Type> interpNgbProcCellData(0);
        if (ptProcCells.found(pointI))
        {
            const List<labelPair>& pc = ptProcCells[pointI];

            interpNgbProcCellData.setSize(pc.size());

            forAll(pc, cI)
            {
                interpNgbProcCellData[cI] =
                    procCellVfI
                    [
                        pc[cI].first()
                    ]
                    [
                        pc[cI].second()
                    ];
            }
        }

        Field<Type> interpNgbProcBndFaceData(0);
        if (ptProcBndFaces[pointI].size())
//         if (ptProcBndFaces.found(pointI))
        {
            const List<labelPair>& pf = ptProcBndFaces[pointI];

            interpNgbProcBndFaceData.setSize(pf.size());

            forAll(pf, fI)
            {
                interpNgbProcBndFaceData[fI] =
                    procBndFaceVf
                    [
                        pf[fI].first()
                    ]
                    [
                        pf[fI].second()
                    ];
            }
        }

        Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
        Field<Type> source
        (
            interpCells.size()
          + interpBndFaces.size()
          + interpCyclicFaces.size()
          + interpProcFaces.size()
          + glInterpNgbProcBndFaceData.size()
          + glInterpNgbProcCellData.size()
          + interpNgbProcCellData.size()
          + interpNgbProcBndFaceData.size(),
            pTraits<Type>::zero
        );

        label pointID = 0;

        Type avg = pTraits<Type>::zero;

        for (label i=0; i<interpCells.size(); i++)
        {
            source[pointID] = vfI[interpCells[i]];
            avg += sqr(W[pointID])*vfI[interpCells[i]];
            pointID++;
        }

        for (label i=0; i<interpBndFaces.size(); i++)
        {
            label faceID = interpBndFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID] = vf.boundaryField()[patchID][localFaceID];
            avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
            pointID++;
        }

        for (label i=0; i<interpCyclicFaces.size(); i++)
        {
            label faceID = interpCyclicFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            const unallocLabelList& faceCells =
                mesh().boundary()[patchID].faceCells();

            label sizeby2 = faceCells.size()/2;

            if (localFaceID < sizeby2)
            {
                source[pointID] = vfI[faceCells[localFaceID + sizeby2]];
                avg += sqr(W[pointID])*vfI[faceCells[localFaceID + sizeby2]];
                pointID++;
            }
            else
            {
                source[pointID] = vfI[faceCells[localFaceID - sizeby2]];
                avg += sqr(W[pointID])*vfI[faceCells[localFaceID - sizeby2]];
                pointID++;
            }
        }

        for (label i=0; i<interpProcFaces.size(); i++)
        {
            label faceID = interpProcFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID] = vf.boundaryField()[patchID][localFaceID];
            avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
            pointID++;
        }

        for (label i=0; i<glInterpNgbProcBndFaceData.size(); i++)
        {
            source[pointID] = glInterpNgbProcBndFaceData[i];
            avg += sqr(W[pointID])*glInterpNgbProcBndFaceData[i];
            pointID++;
        }

        for (label i=0; i<glInterpNgbProcCellData.size(); i++)
        {
            source[pointID] = glInterpNgbProcCellData[i];
            avg += sqr(W[pointID])*glInterpNgbProcCellData[i];
            pointID++;
        }

        for (label i=0; i<interpNgbProcCellData.size(); i++)
        {
            source[pointID] = interpNgbProcCellData[i];
            avg += sqr(W[pointID])*interpNgbProcCellData[i];
            pointID++;
        }

        for (label i=0; i<interpNgbProcBndFaceData.size(); i++)
        {
            source[pointID] = interpNgbProcBndFaceData[i];
            avg += sqr(W[pointID])*interpNgbProcBndFaceData[i];
            pointID++;
        }

        if (mag(mirrorPlaneTransformation()[pointI].first()) > SMALL)
//         if (mirrorPlaneTransformation().found(pointI))
        {
            const tensor& T = mirrorPlaneTransformation()[pointI].second();

            label oldSize = source.size();

            source.setSize(2*oldSize);

            for (label i=oldSize; i<source.size(); i++)
            {
                source[i] = transform(T, source[i-oldSize]);
            }

            avg += transform(T, avg);
        }

//         avg /= source.size() + SMALL;
        avg /= sum(sqr(W));

        source -= avg;

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

        vector dr = (points[pointI] - o[pointI])/L[pointI];

        ppf[pI] =
            avg
          + coeffs[0]*dr.x()
          + coeffs[1]*dr.y()
          + coeffs[2]*dr.z();
    }

    return tppf;
}


// template<class Type>
// tmp<Field<Type> > newLeastSquaresVolPointInterpolation::interpolate
// (
//     const polyPatch& patch,
//     const GeometricField<Type, fvPatchField, volMesh>& vf
// ) const
// {
//     if (debug)
//     {
//         Info<< "newLeastSquaresVolPointInterpolation::interpolate("
//             << "const GeometricField<Type, fvPatchField, volMesh>&, "
//             << "GeometricField<Type, pointPatchField, pointMesh>&) : "
//             << "interpolating field from cells to points"
//             << endl;
//     }

//     Info << "patch cell to point interpolation" << endl;

// #if (defined(OPENFOAM) || defined(OPENFOAMESIORFOUNDATION))
//     const Field<Type>& vfI = vf.primitiveField();
// #else
//     const Field<Type>& vfI = vf.internalField();
// #endif

//     tmp<Field<Type> > tppf
//     (
//         new Field<Type>
//         (
//             patch.nPoints(),
//             pTraits<Type>::zero
//         )
//     );
// #if (defined(OPENFOAM) || defined(OPENFOAMESIORFOUNDATION))
//     Field<Type>& ppf = tppf.ref();
// #else
//     Field<Type>& ppf = tppf();
// #endif

//     const labelList& meshPoints = patch.meshPoints();

//     const vectorField& points = mesh().points();

//     label nCoeffs = 3;
//     const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();
//     const labelListList& ptCells = mesh().pointCells();
//     const labelListList& ptBndFaces = pointBndFaces();
//     const labelListList& ptProcFaces = pointProcFaces();

//     const FieldField<Field, scalar>& w = weights();
//     const vectorField& o = origins();

//     Map<Field<Type> > ptNgbProcBndFaceFieldData;
//     pointNgbProcBndFaceFieldData(vf, ptNgbProcBndFaceFieldData);

//     Map<Field<Type> > gPtNgbProcCellFieldData;
//     globalPointNgbProcCellFieldData(vf, gPtNgbProcCellFieldData);

//     const Map<List<labelPair> >& ptProcCells = pointProcCells();

//     FieldField<Field, Type> procCellVfI = procCellsFieldData(vfI);

//     forAll(ppf, pI)
//     {
//         label pointI = meshPoints[pI];

//         const scalarRectangularMatrix& curInvMatrix = invMatrices[pointI];

//         const scalarField& W = w[pointI];

//         const labelList& interpCells = ptCells[pointI];
//         const labelList& interpBndFaces = ptBndFaces[pointI];
//         const labelList& interpProcFaces = ptProcFaces[pointI];

//         Field<Type> interpNgbProcBndFaceData(0);
//         if (ptNgbProcBndFaceFieldData.found(pointI))
//         {
//             interpNgbProcBndFaceData = ptNgbProcBndFaceFieldData[pointI];
//         }

//         Field<Type> glInterpNgbProcCellData(0);
//         if (gPtNgbProcCellFieldData.found(pointI))
//         {
//             glInterpNgbProcCellData = gPtNgbProcCellFieldData[pointI];
//         }

//         Field<Type> interpNgbProcCellData(0);
//         if (ptProcCells.found(pointI))
//         {
//             const List<labelPair>& pc = ptProcCells[pointI];

//             interpNgbProcCellData.setSize(pc.size());

//             forAll(pc, cI)
//             {
//                 interpNgbProcCellData[cI] =
//                     procCellVfI
//                     [
//                         pc[cI].first()
//                     ]
//                     [
//                         pc[cI].second()
//                     ];
//             }
//         }

//         Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
//         Field<Type> source
//         (
//             interpCells.size()
//           + interpBndFaces.size()
//           + interpProcFaces.size()
//           + interpNgbProcBndFaceData.size()
//           + glInterpNgbProcCellData.size()
//           + interpNgbProcCellData.size(),
//             pTraits<Type>::zero
//         );

//         label pointID = 0;

//         Type avg = pTraits<Type>::zero;

//         for (label i=0; i<interpCells.size(); i++)
//         {
//             source[pointID] = vfI[interpCells[i]];
//             avg += sqr(W[pointID])*vfI[interpCells[i]];
//             pointID++;
//         }

//         for (label i=0; i<interpBndFaces.size(); i++)
//         {
//             label faceID = interpBndFaces[i];
//             label patchID =
//                 mesh().boundaryMesh().whichPatch(faceID);

//             label start = mesh().boundaryMesh()[patchID].start();
//             label localFaceID = faceID - start;

//             source[pointID] = vf.boundaryField()[patchID][localFaceID];
//             avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
//             pointID++;
//         }

//         for (label i=0; i<interpProcFaces.size(); i++)
//         {
//             label faceID = interpProcFaces[i];
//             label patchID =
//                 mesh().boundaryMesh().whichPatch(faceID);

//             label start = mesh().boundaryMesh()[patchID].start();
//             label localFaceID = faceID - start;

//             source[pointID] = vf.boundaryField()[patchID][localFaceID];
//             avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
//             pointID++;
//         }

//         for (label i=0; i<interpNgbProcBndFaceData.size(); i++)
//         {
//             source[pointID] = interpNgbProcBndFaceData[i];
//             avg += sqr(W[pointID])*interpNgbProcBndFaceData[i];
//             pointID++;
//         }

//         for (label i=0; i<glInterpNgbProcCellData.size(); i++)
//         {
//             source[pointID] = glInterpNgbProcCellData[i];
//             avg += sqr(W[pointID])*glInterpNgbProcCellData[i];
//             pointID++;
//         }

//         for (label i=0; i<interpNgbProcCellData.size(); i++)
//         {
//             source[pointID] = interpNgbProcCellData[i];
//             avg += sqr(W[pointID])*interpNgbProcCellData[i];
//             pointID++;
//         }

//         if (mirrorPlaneTransformation().found(pointI))
//         {
//             const tensor& T = mirrorPlaneTransformation()[pointI].second();

//             label oldSize = source.size();

//             source.setSize(2*oldSize);

//             for (label i=oldSize; i<source.size(); i++)
//             {
//                 source[i] = transform(T, source[i-oldSize]);
//             }

//             avg += transform(T, avg);
//         }

// //         avg /= source.size() + SMALL;
//         avg /= sum(sqr(W));

//         source -= avg;

//         for (label i=0; i<nCoeffs; i++)
//         {
//             for (label j=0; j<source.size(); j++)
//             {
//                 coeffs[i] += curInvMatrix[i][j]*source[j];
//             }
//         }

//         vector dr = points[pointI] - o[pointI];

//         ppf[pI] =
//             avg
//           + coeffs[0]*dr.x()
//           + coeffs[1]*dr.y()
//           + coeffs[2]*dr.z();
//     }

//     return tppf;
// }


template<class Type>
Type newLeastSquaresVolPointInterpolation::interpolate
(
    const label pointIndex,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::interpolate("
            << "const label, "
            << "const GeometricField<Type, fvPatchField, volMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    const Field<Type>& vfI = vf.primitiveField();

    Type pf = pTraits<Type>::zero;

    const vectorField& points = mesh().points();

    const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();

    label nCoeffs = 3;
    const labelListList& ptCells = mesh().pointCells();
    const labelListList& ptBndFaces = pointBndFaces();
    const labelListList& ptCyclicFaces = pointCyclicFaces();
    const labelListList& ptProcFaces = pointProcFaces();

    Map<Field<Type> > gPtNgbProcBndFaceFieldData;
    globalPointNgbProcBndFaceFieldData(vf, gPtNgbProcBndFaceFieldData);

    Map<Field<Type> > gPtNgbProcCellFieldData;
    globalPointNgbProcCellFieldData(vf, gPtNgbProcCellFieldData);

    const Map<List<labelPair> >& ptProcCells = pointProcCells();

    const List<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();
//     const Map<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();

    const FieldField<Field, scalar>& w = weights();
    const vectorField& o = origins();
    const scalarField& L = refL();

    FieldField<Field, Type> procCellVfI = procCellsFieldData(vfI);
    FieldField<Field, Type> procBndFaceVf = procBndFacesFieldData(vf);

    // Calc point field value
    {
        label pointI = pointIndex;

        const scalarRectangularMatrix& curInvMatrix = invMatrices[pointI];

        const scalarField& W = w[pointI];

        const labelList& interpCells = ptCells[pointI];
        const labelList& interpBndFaces = ptBndFaces[pointI];
        const labelList& interpCyclicFaces = ptCyclicFaces[pointI];
        const labelList& interpProcFaces = ptProcFaces[pointI];

        Field<Type> glInterpNgbProcBndFaceData(0);
        if (gPtNgbProcBndFaceFieldData.found(pointI))
        {
            glInterpNgbProcBndFaceData = gPtNgbProcBndFaceFieldData[pointI];
        }

        Field<Type> glInterpNgbProcCellData(0);
        if (gPtNgbProcCellFieldData.found(pointI))
        {
            glInterpNgbProcCellData = gPtNgbProcCellFieldData[pointI];
        }

        Field<Type> interpNgbProcCellData(0);
        if (ptProcCells.found(pointI))
        {
            const List<labelPair>& pc = ptProcCells[pointI];

            interpNgbProcCellData.setSize(pc.size());

            forAll(pc, cI)
            {
                interpNgbProcCellData[cI] =
                    procCellVfI
                    [
                        pc[cI].first()
                    ]
                    [
                        pc[cI].second()
                    ];
            }
        }

        Field<Type> interpNgbProcBndFaceData(0);
        if (ptProcBndFaces[pointI].size())
//         if (ptProcBndFaces.found(pointI))
        {
            const List<labelPair>& pf = ptProcBndFaces[pointI];

            interpNgbProcBndFaceData.setSize(pf.size());

            forAll(pf, fI)
            {
                interpNgbProcBndFaceData[fI] =
                    procBndFaceVf
                    [
                        pf[fI].first()
                    ]
                    [
                        pf[fI].second()
                    ];
            }
        }

        Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
        Field<Type> source
        (
            interpCells.size()
          + interpBndFaces.size()
          + interpCyclicFaces.size()
          + interpProcFaces.size()
          + glInterpNgbProcBndFaceData.size()
          + glInterpNgbProcCellData.size()
          + interpNgbProcCellData.size()
          + interpNgbProcBndFaceData.size(),
            pTraits<Type>::zero
        );

        label pointID = 0;

        Type avg = pTraits<Type>::zero;

        for (label i=0; i<interpCells.size(); i++)
        {
            source[pointID] = vfI[interpCells[i]];
            avg += sqr(W[pointID])*vfI[interpCells[i]];
            pointID++;
        }

        for (label i=0; i<interpBndFaces.size(); i++)
        {
            label faceID = interpBndFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID] = vf.boundaryField()[patchID][localFaceID];
            avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
            pointID++;
        }

        for (label i=0; i<interpCyclicFaces.size(); i++)
        {
            label faceID = interpCyclicFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            const unallocLabelList& faceCells =
                mesh().boundary()[patchID].faceCells();

            label sizeby2 = faceCells.size()/2;

            if (localFaceID < sizeby2)
            {
                source[pointID] = vfI[faceCells[localFaceID + sizeby2]];
                avg += sqr(W[pointID])*vfI[faceCells[localFaceID + sizeby2]];
                pointID++;
            }
            else
            {
                source[pointID] = vfI[faceCells[localFaceID - sizeby2]];
                avg += sqr(W[pointID])*vfI[faceCells[localFaceID - sizeby2]];
                pointID++;
            }
        }

        for (label i=0; i<interpProcFaces.size(); i++)
        {
            label faceID = interpProcFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID] = vf.boundaryField()[patchID][localFaceID];
            avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
            pointID++;
        }

        for (label i=0; i<glInterpNgbProcBndFaceData.size(); i++)
        {
            source[pointID] = glInterpNgbProcBndFaceData[i];
            avg += sqr(W[pointID])*glInterpNgbProcBndFaceData[i];
            pointID++;
        }

        for (label i=0; i<glInterpNgbProcCellData.size(); i++)
        {
            source[pointID] = glInterpNgbProcCellData[i];
            avg += sqr(W[pointID])*glInterpNgbProcCellData[i];
            pointID++;
        }

        for (label i=0; i<interpNgbProcCellData.size(); i++)
        {
            source[pointID] = interpNgbProcCellData[i];
            avg += sqr(W[pointID])*interpNgbProcCellData[i];
            pointID++;
        }

        for (label i=0; i<interpNgbProcBndFaceData.size(); i++)
        {
            source[pointID] = interpNgbProcBndFaceData[i];
            avg += sqr(W[pointID])*interpNgbProcBndFaceData[i];
            pointID++;
        }

        if (mag(mirrorPlaneTransformation()[pointI].first()) > SMALL)
//         if (mirrorPlaneTransformation().found(pointI))
        {
            const tensor& T = mirrorPlaneTransformation()[pointI].second();

            label oldSize = source.size();

            source.setSize(2*oldSize);

            for (label i=oldSize; i<source.size(); i++)
            {
                source[i] = transform(T, source[i-oldSize]);
            }

            avg += transform(T, avg);
        }

//         avg /= source.size() + SMALL;
        avg /= sum(sqr(W));

        source -= avg;

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

        vector dr = (points[pointI] - o[pointI])/L[pointI];

        pf =
            avg
          + coeffs[0]*dr.x()
          + coeffs[1]*dr.y()
          + coeffs[2]*dr.z();
    }

    return pf;
}


// template<class Type>
// Type newLeastSquaresVolPointInterpolation::interpolate
// (
//     const label pointIndex,
//     const GeometricField<Type, fvPatchField, volMesh>& vf
// ) const
// {
//     if (debug)
//     {
//         Info<< "newLeastSquaresVolPointInterpolation::interpolate("
//             << "const label, "
//             << "const GeometricField<Type, fvPatchField, volMesh>&) : "
//             << "interpolating field from cells to points"
//             << endl;
//     }

// #if (defined(OPENFOAM) || defined(OPENFOAMESIORFOUNDATION))
//     const Field<Type>& vfI = vf.primitiveField();
// #else
//     const Field<Type>& vfI = vf.internalField();
// #endif
//     Type pf = pTraits<Type>::zero;


//     const vectorField& points = mesh().points();

//     label nCoeffs = 3;
//     const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();
//     const labelListList& ptCells = mesh().pointCells();
//     const labelListList& ptBndFaces = pointBndFaces();
//     const labelListList& ptProcFaces = pointProcFaces();

//     const FieldField<Field, scalar>& w = weights();
//     const vectorField& o = origins();

//     Map<Field<Type> > ptNgbProcBndFaceFieldData;
//     pointNgbProcBndFaceFieldData(vf, ptNgbProcBndFaceFieldData);

//     Map<Field<Type> > gPtNgbProcCellFieldData;
//     globalPointNgbProcCellFieldData(vf, gPtNgbProcCellFieldData);

//     const Map<List<labelPair> >& ptProcCells = pointProcCells();

//     FieldField<Field, Type> procCellVfI = procCellsFieldData(vfI);

//     // Calc point field value
//     {
//         label pointI = pointIndex;

//         const scalarRectangularMatrix& curInvMatrix = invMatrices[pointI];

//         const scalarField& W = w[pointI];

//         const labelList& interpCells = ptCells[pointI];
//         const labelList& interpBndFaces = ptBndFaces[pointI];
//         const labelList& interpProcFaces = ptProcFaces[pointI];

//         Field<Type> interpNgbProcBndFaceData(0);
//         if (ptNgbProcBndFaceFieldData.found(pointI))
//         {
//             interpNgbProcBndFaceData = ptNgbProcBndFaceFieldData[pointI];
//         }

//         Field<Type> glInterpNgbProcCellData(0);
//         if (gPtNgbProcCellFieldData.found(pointI))
//         {
//             glInterpNgbProcCellData = gPtNgbProcCellFieldData[pointI];
//         }

//         Field<Type> interpNgbProcCellData(0);
//         if (ptProcCells.found(pointI))
//         {
//             const List<labelPair>& pc = ptProcCells[pointI];

//             interpNgbProcCellData.setSize(pc.size());

//             forAll(pc, cI)
//             {
//                 interpNgbProcCellData[cI] =
//                     procCellVfI
//                     [
//                         pc[cI].first()
//                     ]
//                     [
//                         pc[cI].second()
//                     ];
//             }
//         }

//         Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
//         Field<Type> source
//         (
//             interpCells.size()
//           + interpBndFaces.size()
//           + interpProcFaces.size()
//           + interpNgbProcBndFaceData.size()
//           + glInterpNgbProcCellData.size()
//           + interpNgbProcCellData.size(),
//             pTraits<Type>::zero
//         );

//         label pointID = 0;

//         Type avg = pTraits<Type>::zero;

//         for (label i=0; i<interpCells.size(); i++)
//         {
//             source[pointID] = vfI[interpCells[i]];
//             avg += sqr(W[pointID])*vfI[interpCells[i]];
//             pointID++;
//         }

//         for (label i=0; i<interpBndFaces.size(); i++)
//         {
//             label faceID = interpBndFaces[i];
//             label patchID =
//                 mesh().boundaryMesh().whichPatch(faceID);

//             label start = mesh().boundaryMesh()[patchID].start();
//             label localFaceID = faceID - start;

//             source[pointID] = vf.boundaryField()[patchID][localFaceID];
//             avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
//             pointID++;
//         }

//         for (label i=0; i<interpProcFaces.size(); i++)
//         {
//             label faceID = interpProcFaces[i];
//             label patchID =
//                 mesh().boundaryMesh().whichPatch(faceID);

//             label start = mesh().boundaryMesh()[patchID].start();
//             label localFaceID = faceID - start;

//             source[pointID] = vf.boundaryField()[patchID][localFaceID];
//             avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
//             pointID++;
//         }

//         for (label i=0; i<interpNgbProcBndFaceData.size(); i++)
//         {
//             source[pointID] = interpNgbProcBndFaceData[i];
//             avg += sqr(W[pointID])*interpNgbProcBndFaceData[i];
//             pointID++;
//         }

//         for (label i=0; i<glInterpNgbProcCellData.size(); i++)
//         {
//             source[pointID] = glInterpNgbProcCellData[i];
//             avg += sqr(W[pointID])*glInterpNgbProcCellData[i];
//             pointID++;
//         }

//         for (label i=0; i<interpNgbProcCellData.size(); i++)
//         {
//             source[pointID] = interpNgbProcCellData[i];
//             avg += sqr(W[pointID])*interpNgbProcCellData[i];
//             pointID++;
//         }

//         if (mirrorPlaneTransformation().found(pointI))
//         {
//             const tensor& T = mirrorPlaneTransformation()[pointI].second();

//             label oldSize = source.size();

//             source.setSize(2*oldSize);

//             for (label i=oldSize; i<source.size(); i++)
//             {
//                 source[i] = transform(T, source[i-oldSize]);
//             }

//             avg += transform(T, avg);
//         }

// //         avg /= source.size() + SMALL;
//         avg /= sum(sqr(W));

//         source -= avg;

//         for (label i=0; i<nCoeffs; i++)
//         {
//             for (label j=0; j<source.size(); j++)
//             {
//                 coeffs[i] += curInvMatrix[i][j]*source[j];
//             }
//         }

//         vector dr = points[pointI] - o[pointI];

//         pf =
//             avg
//           + coeffs[0]*dr.x()
//           + coeffs[1]*dr.y()
//           + coeffs[2]*dr.z();
//     }

//     return pf;
// }


// template<class Type>
// void newLeastSquaresVolPointInterpolation::pointNgbProcBndFaceFieldData
// (
//     const GeometricField<Type, fvPatchField, volMesh>& vf,
//     Map<Field<Type> >& fieldData
// ) const
// {
//     if (debug)
//     {
//         Info<< "newLeastSquaresVolPointInterpolation::"
//             << "pointNgbProcBndFaceFieldData("
//             << "const GeometricField<Type, fvPatchField, volMesh>&) : "
//             << "extracting bnd face data from ngb processors"
//             << endl;
//     }

//     const labelListList& ptBndFaces = pointBndFaces();

//     if (Pstream::parRun())
//     {
//         forAll(mesh().boundaryMesh(), patchI)
//         {
//             if
//             (
//                 mesh().boundaryMesh()[patchI].type()
//              == processorPolyPatch::typeName
//             )
//             {
//                 const processorPolyPatch& procPatch =
//                     refCast<const processorPolyPatch>
//                     (
//                         mesh().boundaryMesh()[patchI]
//                     );

//                 const labelList& bndPoints = procPatch.boundaryPoints();
//                 const labelList& mshPoints = procPatch.meshPoints();

//                 FieldField<Field, Type> bndPointsFaceData
//                 (
//                     bndPoints.size()
//                 );

//                 forAll(bndPoints, pointI)
//                 {
//                     label curPoint = bndPoints[pointI];
//                     label curMeshPoint = mshPoints[curPoint];

//                     const labelList& curPointBndFaces =
//                         ptBndFaces[curMeshPoint];

//                     Field<Type> curPointFaceData(curPointBndFaces.size());

//                     forAll(curPointFaceData, faceI)
//                     {
//                         label faceID = curPointBndFaces[faceI];
//                         label patchID =
//                             mesh().boundaryMesh().whichPatch(faceID);

//                         label start = mesh().boundaryMesh()[patchID].start();
//                         label localFaceID = faceID - start;

//                         curPointFaceData[faceI] =
//                             vf.boundaryField()[patchID][localFaceID];
//                     }

//                     bndPointsFaceData.set
//                     (
//                         pointI,
//                         new Field<Type>(curPointFaceData)
//                     );
//                 }

//                 // Parallel data exchange
//                 {
//                     OPstream toNeighbProc
//                     (
//                         Pstream::blocking,
//                         procPatch.neighbProcNo()
//                         // size of field
//                     );

//                     toNeighbProc << bndPoints << bndPointsFaceData;
//                 }

//                 FieldField<Field, Type> ngbBndPointsFaceData
//                 (
//                     bndPoints.size()
//                 );

//                 labelList ngbBndPoints(bndPoints.size());

//                 {
//                     IPstream fromNeighbProc
//                     (
//                         Pstream::blocking,
//                         procPatch.neighbProcNo()
//                         // size of field
//                     );

//                     fromNeighbProc >> ngbBndPoints >> ngbBndPointsFaceData;
//                 }

//                 const labelList& glPoints =
//                     mesh().globalData().sharedPointLabels();

//                 forAll(bndPoints, pointI)
//                 {
//                     label curPoint = bndPoints[pointI];
//                     label curMeshPoint = mshPoints[curPoint];

//                     label gpIndex = findIndex(glPoints, curMeshPoint);

//                     // non-global points
//                     if (gpIndex == -1)
//                     {
//                         label curNgbPoint = procPatch.neighbPoints()[curPoint];

//                         label curNgbBndPoint =
//                             findIndex(ngbBndPoints, curNgbPoint);

//                         fieldData.insert
//                         (
//                             curMeshPoint,
//                             ngbBndPointsFaceData[curNgbBndPoint]
//                         );
//                     }
//                 }
//             }
//         }

//         // Global boundary points
//         if (mesh().globalData().nGlobalPoints())
//         {
//             const labelList& spLabels =
//                 mesh().globalData().sharedPointLabels();

//             const labelList& spAddressing =
//                 mesh().globalData().sharedPointAddr();

//             for (label k=0; k<mesh().globalData().nGlobalPoints(); k++)
//             {
//                 List<List<Type> > procBndFaceData(Pstream::nProcs());

//                 label curSpIndex = findIndex(spAddressing, k);

//                 if (curSpIndex != -1)
//                 {
//                     label curMeshPoint = spLabels[curSpIndex];

//                     const labelList& curBndFaces = ptBndFaces[curMeshPoint];

//                     procBndFaceData[Pstream::myProcNo()] =
//                         List<Type>(curBndFaces.size());

//                     forAll (curBndFaces, faceI)
//                     {
//                         label faceID = curBndFaces[faceI];
//                         label patchID =
//                             mesh().boundaryMesh().whichPatch(faceID);

//                         label start = mesh().boundaryMesh()[patchID].start();
//                         label localFaceID = faceID - start;

//                         procBndFaceData[Pstream::myProcNo()][faceI] =
//                             vf.boundaryField()[patchID][localFaceID];
//                     }
//                 }
//                 else
//                 {
//                     procBndFaceData[Pstream::myProcNo()] = List<Type>(0);
//                 }

//                 Pstream::gatherList(procBndFaceData);
//                 Pstream::scatterList(procBndFaceData);

//                 if (curSpIndex != -1)
//                 {
//                     label curMeshPoint = spLabels[curSpIndex];

//                     label nAllFaces = 0;
//                     forAll(procBndFaceData, procI)
//                     {
//                         if (procI != Pstream::myProcNo())
//                         {
//                             nAllFaces += procBndFaceData[procI].size();
//                         }
//                     }

//                     Field<Type> allFaceData(nAllFaces, pTraits<Type>::zero);

//                     label counter = 0;
//                     forAll(procBndFaceData, procI)
//                     {
//                         if (procI != Pstream::myProcNo())
//                         {
//                             forAll(procBndFaceData[procI], faceI)
//                             {
//                                 allFaceData[counter++] =
//                                     procBndFaceData[procI][faceI];
//                             }
//                         }
//                     }

//                     fieldData.insert
//                     (
//                         curMeshPoint,
//                         allFaceData
//                     );
//                 }
//             }
//         }
//     }
// }

template<class Type>
void newLeastSquaresVolPointInterpolation::globalPointNgbProcBndFaceFieldData
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    Map<Field<Type> >& fieldData
) const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::"
            << "globalPointNgbProcBndFaceFieldData("
            << "const GeometricField<Type, fvPatchField, volMesh>&) : "
            << "extracting bnd face data from ngb processors"
            << endl;
    }

//     const labelListList& ptCells = mesh().pointCells();

//     const Field<Type>& vfI = vf.internalField();

    if (Pstream::parRun())
    {
        //notImplemented("globalPointNgbProcBndFaceFieldData");
        if (processorBoundariesExist_)
        {
            FatalErrorIn
            (
                "void newLeastSquaresVolPointInterpolation::"
                "globalPointNgbProcBndFaceFieldData"
            )   << "Currently only the 'none' decomposition is allowed in "
                << "parallel" << abort(FatalError);
        }
    }
}

template<class Type>
void newLeastSquaresVolPointInterpolation::globalPointNgbProcCellFieldData
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    Map<Field<Type> >& fieldData
) const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::"
            << "globalPointNgbProcCellFieldData("
            << "const GeometricField<Type, fvPatchField, volMesh>&) : "
            << "extracting cell data from ngb processors"
            << endl;
    }

    if (Pstream::parRun())
    {
        // notImplemented("globalPointNgbProcCellFieldData");
        if (processorBoundariesExist_)
        {
            FatalErrorIn
            (
                "void newLeastSquaresVolPointInterpolation::"
                "globalPointNgbProcCellFieldData"
            )   << "Currently only the 'none' decomposition is allowed in "
                << "parallel" << abort(FatalError);
        }
    }
}


template<class Type>
Foam::tmp<Foam::FieldField<Foam::Field, Type> >
newLeastSquaresVolPointInterpolation::procCellsFieldData
(
    const Field<Type>& psi
) const
{
    tmp<FieldField<Field, Type> > tprocPsi
    (
        new FieldField<Field, Type>(Pstream::nProcs())
    );
    FieldField<Field, Type>& procPsi = tprocPsi.ref();

    forAll (procPsi, procI)
    {
        procPsi.set
        (
            procI,
            new Field<Type>
            (
                procCellCentres()[procI].size(),
                pTraits<Type>::zero
            )
        );
    }

    if (Pstream::parRun())
    {
        //notImplemented("procCellsFieldData");
        if (processorBoundariesExist_)
        {
            FatalErrorIn
            (
                "void newLeastSquaresVolPointInterpolation::"
                "procCellsFieldData"
            )   << "Currently only the 'none' decomposition is allowed in "
                << "parallel" << abort(FatalError);
        }
    }

    return tprocPsi;
}


template<class Type>
Foam::tmp<Foam::FieldField<Foam::Field, Type> >
newLeastSquaresVolPointInterpolation::procBndFacesFieldData
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<FieldField<Field, Type> > tprocPsi
    (
        new FieldField<Field, Type>(Pstream::nProcs())
    );
    FieldField<Field, Type>& procPsi = tprocPsi.ref();

    forAll (procPsi, procI)
    {
        procPsi.set
        (
            procI,
            new Field<Type>
            (
                procBndFaceCentres()[procI].size(),
                pTraits<Type>::zero
            )
        );
    }

    if (Pstream::parRun())
    {
        // notImplemented("procBndFacesFieldData");
        if (processorBoundariesExist_)
        {
            FatalErrorIn
            (
                "void newLeastSquaresVolPointInterpolation::"
                "procBndFacesFieldData"
            )   << "Currently only the 'none' decomposition is allowed in "
                << "parallel" << abort(FatalError);
        }
    }

    return tprocPsi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

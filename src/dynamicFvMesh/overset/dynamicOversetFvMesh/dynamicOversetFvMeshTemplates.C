/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "volFields.H"
#include "fvMatrix.H"
#include "cellCellStencilObject.H"
#include "oversetFvPatchField.H"
#include "calculatedProcessorFvPatchField.H"
#include "lduInterfaceFieldPtrsList.H"
#include "processorFvPatch.H"
#include "syncTools.H"

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::dynamicOversetFvMesh::interpolate(Field<T>& psi) const
{
    const cellCellStencil& overlap = Stencil::New(*this);
    const labelListList& stencil = overlap.cellStencil();

    if (stencil.size() != nCells())
    {
        return;
    }

    const mapDistribute& map = overlap.cellInterpolationMap();
    const List<scalarList>& wghts = overlap.cellInterpolationWeights();
    const labelList& cellIDs = overlap.interpolationCells();
    const scalarList& factor = overlap.cellInterpolationWeight();

    Field<T> work(psi);
    map.mapDistributeBase::distribute(work, UPstream::msgType()+1);

    forAll(cellIDs, i)
    {
        label celli = cellIDs[i];

        const scalarList& w = wghts[celli];
        const labelList& nbrs = stencil[celli];
        const scalar f = factor[celli];

        T s(pTraits<T>::zero);
        forAll(nbrs, nbrI)
        {
            s += w[nbrI]*work[nbrs[nbrI]];
        }
        //Pout<< "Interpolated value:" << s << endl;
        //T oldPsi = psi[celli];
        psi[celli] = (1.0-f)*psi[celli] + f*s;
        //Pout<< "psi was:" << oldPsi << " now:" << psi[celli] << endl;
    }
}


template<class GeoField>
void Foam::dynamicOversetFvMesh::interpolate(GeoField& psi) const
{
    interpolate(psi.primitiveFieldRef());
    psi.correctBoundaryConditions();
}


template<class GeoField>
void Foam::dynamicOversetFvMesh::interpolate(const wordHashSet& suppressed)
{
    HashTable<GeoField*> flds(this->objectRegistry::lookupClass<GeoField>());
    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        const word& name = iter()->name();
        if (!suppressed.found(baseName(name)))
        {
            if (debug)
            {
                Pout<< "dynamicOversetFvMesh::interpolate: interpolating : "
                    << name << endl;
            }
            interpolate(iter()->primitiveFieldRef());
        }
        else
        {
            if (debug)
            {
                Pout<< "dynamicOversetFvMesh::interpolate: skipping : " << name
                    << endl;
            }
        }
    }
}


template<class GeoField, class PatchType>
void Foam::dynamicOversetFvMesh::correctBoundaryConditions
(
    typename GeoField::Boundary& bfld,
    const bool typeOnly
)
{
    const label nReq = Pstream::nRequests();

    forAll(bfld, patchi)
    {
        if (typeOnly == isA<PatchType>(bfld[patchi]))
        {
            bfld[patchi].initEvaluate(Pstream::defaultCommsType);
        }
    }

    // Block for any outstanding requests
    if (Pstream::parRun())
    {
        Pstream::waitRequests(nReq);
    }

    forAll(bfld, patchi)
    {
        if (typeOnly == isA<PatchType>(bfld[patchi]))
        {
            bfld[patchi].evaluate(Pstream::defaultCommsType);
        }
    }
}


template<class Type>
Foam::tmp<Foam::scalarField> Foam::dynamicOversetFvMesh::normalisation
(
    const fvMatrix<Type>& m
) const
{
    // Determine normalisation. This is normally the original diagonal plus
    // remote contributions. This needs to be stabilised for hole cells
    // which can have a zero diagonal. Assume that if any component has
    // a non-zero diagonal the cell does not need stabilisation.
    tmp<scalarField> tnorm(new scalarField(m.diag()));
    scalarField& norm = tnorm.ref();

    // Add remote coeffs to duplicate behaviour of fvMatrix::addBoundaryDiag
    fvMatrix<Type>& mRef = const_cast<fvMatrix<Type>&>(m);
    const FieldField<Field, Type>& internalCoeffs = mRef.internalCoeffs();
    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        forAll(internalCoeffs, patchi)
        {
            const labelUList& fc = lduAddr().patchAddr(patchi);
            const Field<Type>& intCoeffs = internalCoeffs[patchi];
            const scalarField cmptCoeffs(intCoeffs.component(cmpt));
            forAll(fc, i)
            {
                norm[fc[i]] += cmptCoeffs[i];
            }
        }
    }

    // Count number of problematic cells
    label nZeroDiag = 0;
    forAll(norm, celli)
    {
        const scalar& n = norm[celli];
        if (magSqr(n) < sqr(SMALL))
        {
            //Pout<< "For field " << m.psi().name()
            //    << " have diagonal " << n << " for cell " << celli
            //    << " at:" << cellCentres()[celli] << endl;
            nZeroDiag++;
        }
    }

    reduce(nZeroDiag, sumOp<label>());

    if (debug)
    {
        Pout<< "For field " << m.psi().name() << " have zero diagonals for "
            << nZeroDiag << " cells" << endl;
    }

    if (nZeroDiag > 0)
    {
        // Walk out the norm across hole cells

        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();
        const cellCellStencilObject& overlap = Stencil::New(*this);
        const labelUList& types = overlap.cellTypes();

        label nHoles = 0;
        scalarField extrapolatedNorm(norm);
        forAll(types, celli)
        {
            if (types[celli] == cellCellStencil::HOLE)
            {
                extrapolatedNorm[celli] = -GREAT;
                nHoles++;
            }
        }

        PackedBoolList isFront(nFaces());
        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            label ownType = types[own[facei]];
            label neiType = types[nei[facei]];
            if
            (
                (ownType == cellCellStencil::HOLE)
             != (neiType == cellCellStencil::HOLE)
            )
            {
                isFront.set(facei);
            }
        }
        labelList nbrTypes;
        syncTools::swapBoundaryCellList(*this, types, nbrTypes);
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            label ownType = types[own[facei]];
            label neiType = nbrTypes[facei-nInternalFaces()];
            if
            (
                (ownType == cellCellStencil::HOLE)
             != (neiType == cellCellStencil::HOLE)
            )
            {
                isFront.set(facei);
            }
        }


        while (true)
        {
            scalarField nbrNorm;
            syncTools::swapBoundaryCellList(*this, extrapolatedNorm, nbrNorm);

            PackedBoolList newIsFront(nFaces());
            scalarField newNorm(extrapolatedNorm);

            label nChanged = 0;
            for (const label facei : isFront)
            {
                if (extrapolatedNorm[own[facei]] == -GREAT)
                {
                    // Average owner cell, add faces to newFront
                    newNorm[own[facei]] = cellAverage
                    (
                        types,
                        nbrTypes,
                        extrapolatedNorm,
                        nbrNorm,
                        own[facei],
                        newIsFront
                    );
                    nChanged++;
                }
                if
                (
                    isInternalFace(facei)
                 && extrapolatedNorm[nei[facei]] == -GREAT
                )
                {
                    // Average nei cell, add faces to newFront
                    newNorm[nei[facei]] = cellAverage
                    (
                        types,
                        nbrTypes,
                        extrapolatedNorm,
                        nbrNorm,
                        nei[facei],
                        newIsFront
                    );
                    nChanged++;
                }
            }

            reduce(nChanged, sumOp<label>());
            if (nChanged == 0)
            {
                break;
            }

            // Transfer new front
            extrapolatedNorm.transfer(newNorm);
            isFront.transfer(newIsFront);
            syncTools::syncFaceList(*this, isFront, maxEqOp<unsigned int>());
        }


        forAll(norm, celli)
        {
            scalar& n = norm[celli];
            if (magSqr(n) < sqr(SMALL))
            {
                //Pout<< "For field " << m.psi().name()
                //    << " for cell " << celli
                //    << " at:" << cellCentres()[celli]
                //    << " have norm " << n
                //    << " have extrapolated norm " << extrapolatedNorm[celli]
                //    << endl;
                // Override the norm
                n = extrapolatedNorm[celli];
            }
        }
    }
    return tnorm;
}


template<class Type>
void Foam::dynamicOversetFvMesh::addInterpolation
(
    fvMatrix<Type>& m,
    const scalarField& normalisation
) const
{
    const cellCellStencilObject& overlap = Stencil::New(*this);
    const List<scalarList>& wghts = overlap.cellInterpolationWeights();
    const labelListList& stencil = overlap.cellStencil();
    const labelList& cellIDs = overlap.interpolationCells();
    const scalarList& factor = overlap.cellInterpolationWeight();
    const labelUList& types = overlap.cellTypes();


    // Force asymmetric matrix (if it wasn't already)
    scalarField& lower = m.lower();
    scalarField& upper = m.upper();
    Field<Type>& source = m.source();
    scalarField& diag = m.diag();


    // Get the addressing. Note that the addressing is now extended with
    // any interpolation faces.
    const lduAddressing& addr = lduAddr();
    const labelUList& upperAddr = addr.upperAddr();
    const labelUList& lowerAddr = addr.lowerAddr();
    const labelUList& ownerStartAddr = addr.ownerStartAddr();
    const labelUList& losortAddr = addr.losortAddr();
    const lduInterfacePtrsList& interfaces = allInterfaces_;

    if (!isA<fvMeshPrimitiveLduAddressing>(addr))
    {
        FatalErrorInFunction
            << "Problem : addressing is not fvMeshPrimitiveLduAddressing"
            << exit(FatalError);
    }



    // 1. Adapt lduMatrix for additional faces and new ordering
    upper.setSize(upperAddr.size(), 0.0);
    inplaceReorder(reverseFaceMap_, upper);
    lower.setSize(lowerAddr.size(), 0.0);
    inplaceReorder(reverseFaceMap_, lower);


    //const label nOldInterfaces = dynamicMotionSolverFvMesh::interfaces().size();
    const label nOldInterfaces = fvMesh::interfaces().size();


    if (interfaces.size() > nOldInterfaces)
    {
        // Extend matrix coefficients
        m.internalCoeffs().setSize(interfaces.size());
        m.boundaryCoeffs().setSize(interfaces.size());

        // 1b. Adapt for additional interfaces
        for
        (
            label patchi = nOldInterfaces;
            patchi < interfaces.size();
            patchi++
        )
        {
            const labelUList& fc = interfaces[patchi].faceCells();
            m.internalCoeffs().set(patchi, new Field<Type>(fc.size(), Zero));
            m.boundaryCoeffs().set(patchi, new Field<Type>(fc.size(), Zero));
        }

        // 1c. Adapt field for additional interfaceFields (note: solver uses
        //     GeometricField::scalarInterfaces() to get hold of interfaces)
        typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

        typename GeoField::Boundary& bfld =
            const_cast<GeoField&>(m.psi()).boundaryFieldRef();

        bfld.setSize(interfaces.size());


        // This gets quite interesting: we do not want to add additional
        // fvPatches (since direct correspondence to polyMesh) so instead
        // add a reference to an existing processor patch
        label addPatchi = 0;
        for (label patchi = 0; patchi < nOldInterfaces; patchi++)
        {
            if (isA<processorFvPatch>(bfld[patchi].patch()))
            {
                addPatchi = patchi;
                break;
            }
        }

        for
        (
            label patchi = nOldInterfaces;
            patchi < interfaces.size();
            patchi++
        )
        {
            bfld.set
            (
                patchi,
                new calculatedProcessorFvPatchField<Type>
                (
                    interfaces[patchi],
                    bfld[addPatchi].patch(),    // dummy processorFvPatch
                    m.psi()
                )
            );
        }
    }


    // 2. Adapt fvMatrix level: faceFluxCorrectionPtr
    // Question: do we need to do this?
    // This seems to be set/used only by the gaussLaplacianScheme and
    // fvMatrix:correction, both of which are outside the linear solver.



    // Clear out existing connections on cells to be interpolated
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: could avoid doing the zeroing of the new faces since these
    //       are set to zero anyway.

    forAll(upperAddr, facei)
    {
        if (types[upperAddr[facei]] == cellCellStencil::INTERPOLATED)
        {
            // Disconnect upper from lower
            label celli = upperAddr[facei];
            lower[facei] *= 1.0-factor[celli];
        }
        if (types[lowerAddr[facei]] == cellCellStencil::INTERPOLATED)
        {
            // Disconnect lower from upper
            label celli = lowerAddr[facei];
            upper[facei] *= 1.0-factor[celli];
        }
    }

    for (label patchi = 0; patchi < nOldInterfaces; ++patchi)
    {
        const labelUList& fc = addr.patchAddr(patchi);
        Field<Type>& intCoeffs = m.internalCoeffs()[patchi];
        Field<Type>& bouCoeffs = m.boundaryCoeffs()[patchi];
        forAll(fc, i)
        {
            label celli = fc[i];
            {
                if (types[celli] == cellCellStencil::INTERPOLATED)
                {
                    scalar f = factor[celli];
                    intCoeffs[i] *= 1.0-f;
                    bouCoeffs[i] *= 1.0-f;
                }
                else if (types[celli] == cellCellStencil::HOLE)
                {
                    intCoeffs[i] = pTraits<Type>::zero;
                    bouCoeffs[i] = pTraits<Type>::zero;
                }
            }
        }
    }



    // Modify matrix
    // ~~~~~~~~~~~~~

    // Do hole cells. Note: maybe put into interpolationCells() loop above?
    forAll(types, celli)
    {
        if (types[celli] == cellCellStencil::HOLE)
        {
            label startLabel = ownerStartAddr[celli];
            label endLabel = ownerStartAddr[celli + 1];

            for (label facei = startLabel; facei < endLabel; facei++)
            {
                upper[facei] = 0.0;
            }

            startLabel = addr.losortStartAddr()[celli];
            endLabel = addr.losortStartAddr()[celli + 1];

            for (label i = startLabel; i < endLabel; i++)
            {
                label facei = losortAddr[i];
                lower[facei] = 0.0;
            }

            diag[celli] = normalisation[celli];
            source[celli] = normalisation[celli]*m.psi()[celli];
        }
    }


    //const globalIndex globalNumbering(V().size());
    //labelList globalCellIDs(overlap.cellInterpolationMap().constructSize());
    //forAll(V(), cellI)
    //{
    //    globalCellIDs[cellI] = globalNumbering.toGlobal(cellI);
    //}
    //overlap.cellInterpolationMap().distribute(globalCellIDs);


    forAll(cellIDs, i)
    {
        label celli = cellIDs[i];

        const scalar f = factor[celli];
        const scalarList& w = wghts[celli];
        const labelList& nbrs = stencil[celli];
        const labelList& nbrFaces = stencilFaces_[celli];
        const labelList& nbrPatches = stencilPatches_[celli];

        if (types[celli] == cellCellStencil::HOLE)
        {
            FatalErrorInFunction << "Found HOLE cell " << celli
                << " at:" << C()[celli]
                << " . Should this be in interpolationCells()????"
                << abort(FatalError);
        }
        else
        {
            // Create interpolation stencil

            diag[celli] *= (1.0-f);
            source[celli] *= (1.0-f);

            forAll(nbrs, nbri)
            {
                label patchi = nbrPatches[nbri];
                label facei = nbrFaces[nbri];

                if (patchi == -1)
                {
                    label nbrCelli = nbrs[nbri];

                    // Add the coefficients
                    const scalar s = normalisation[celli]*f*w[nbri];

                    scalar& u = upper[facei];
                    scalar& l = lower[facei];
                    if (celli < nbrCelli)
                    {
                        diag[celli] += s;
                        u += -s;
                    }
                    else
                    {
                        diag[celli] += s;
                        l += -s;
                    }
                }
                else
                {
                    // Patch face. Store in boundaryCoeffs. Note sign change.
                    //const label globalCelli = globalCellIDs[nbrs[nbri]];
                    //const label proci =
                    //    globalNumbering.whichProcID(globalCelli);
                    //const label remoteCelli =
                    //    globalNumbering.toLocal(proci, globalCelli);
                    //
                    //Pout<< "for cell:" << celli
                    //    << " need weight from remote slot:" << nbrs[nbri]
                    //    << " proc:" << proci << " remote cell:" << remoteCelli
                    //    << " patch:" << patchi
                    //    << " patchFace:" << facei
                    //    << " weight:" << w[nbri]
                    //    << endl;

                    const scalar s = normalisation[celli]*f*w[nbri];
                    m.boundaryCoeffs()[patchi][facei] += pTraits<Type>::one*s;
                    m.internalCoeffs()[patchi][facei] += pTraits<Type>::one*s;

                    // Note: do NOT add to diagonal - this is in the
                    //       internalCoeffs and gets added to the diagonal
                    //       inside fvMatrix::solve
                }
            }

            //if (mag(diag[celli]) < SMALL)
            //{
            //    Pout<< "for cell:" << celli
            //        << " at:" << this->C()[celli]
            //        << " diag:" << diag[celli] << endl;
            //
            //    forAll(nbrs, nbri)
            //    {
            //        label patchi = nbrPatches[nbri];
            //        label facei = nbrFaces[nbri];
            //
            //        const label globalCelli = globalCellIDs[nbrs[nbri]];
            //        const label proci =
            //            globalNumbering.whichProcID(globalCelli);
            //        const label remoteCelli =
            //            globalNumbering.toLocal(proci, globalCelli);
            //
            //        Pout<< " need weight from slot:" << nbrs[nbri]
            //            << " proc:" << proci << " remote cell:"
            //            << remoteCelli
            //            << " patch:" << patchi
            //            << " patchFace:" << facei
            //            << " weight:" << w[nbri]
            //            << endl;
            //    }
            //    Pout<< endl;
            //}
        }
    }
}


template<class Type>
Foam::SolverPerformance<Type> Foam::dynamicOversetFvMesh::solve
(
    fvMatrix<Type>& m,
    const dictionary& dict
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;
    // Check we're running with bcs that can handle implicit matrix manipulation
    typename GeoField::Boundary& bpsi =
        const_cast<GeoField&>(m.psi()).boundaryFieldRef();

    bool hasOverset = false;
    forAll(bpsi, patchi)
    {
        if (isA<oversetFvPatchField<Type>>(bpsi[patchi]))
        {
            hasOverset = true;
            break;
        }
    }

    if (!hasOverset)
    {
        if (debug)
        {
            Pout<< "dynamicOversetFvMesh::solve() :"
                << " bypassing matrix adjustment for field " << m.psi().name()
                << endl;
        }
        //return dynamicMotionSolverFvMesh::solve(m, dict);
        return m.solve(dict);
    }

    if (debug)
    {
        Pout<< "dynamicOversetFvMesh::solve() :"
            << " adjusting matrix for interpolation for field "
            << m.psi().name() << endl;
    }

    // Calculate stabilised diagonal as normalisation for interpolation
    const scalarField norm(normalisation(m));

    if (debug)
    {
        volScalarField scale
        (
            IOobject
            (
                m.psi().name() + "_normalisation",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimensionedScalar(dimless, Zero)
        );
        scale.ref().field() = norm;
        correctBoundaryConditions
        <
            volScalarField,
            oversetFvPatchField<scalar>
        >(scale.boundaryFieldRef(), false);
        scale.write();

        if (debug)
        {
            Pout<< "dynamicOversetFvMesh::solve() :"
                << " writing matrix normalisation for field " << m.psi().name()
                << " to " << scale.name() << endl;
        }
    }


    // Switch to extended addressing (requires mesh::update() having been
    // called)
    active(true);

    // Adapt matrix
    scalarField oldUpper(m.upper());
    scalarField oldLower(m.lower());
    FieldField<Field, Type> oldInt(m.internalCoeffs());
    FieldField<Field, Type> oldBou(m.boundaryCoeffs());
    const label oldSize = bpsi.size();

    addInterpolation(m, norm);

    // Swap psi values so added patches have patchNeighbourField
    correctBoundaryConditions<GeoField, calculatedProcessorFvPatchField<Type>>
    (
        bpsi,
        true
    );


    // Print a bit
    //write(Pout, m, lduAddr(), interfaces());
    //{
    //   const fvSolution& sol = static_cast<const fvSolution&>(*this);
    //    const dictionary& pDict = sol.subDict("solvers").subDict("p");
    //    writeAgglomeration(GAMGAgglomeration::New(m, pDict));
    //}

    // Use lower level solver
    //SolverPerformance<Type> s(dynamicMotionSolverFvMesh::solve(m, dict));
    SolverPerformance<Type> s(m.solve(dict));

    // Restore boundary field
    bpsi.setSize(oldSize);

    // Restore matrix
    m.upper().transfer(oldUpper);
    m.lower().transfer(oldLower);
    m.internalCoeffs().transfer(oldInt);
    m.boundaryCoeffs().transfer(oldBou);

    // Switch to original addressing
    active(false);

    return s;
}


template<class Type>
void Foam::dynamicOversetFvMesh::write
(
    Ostream& os,
    const fvMatrix<Type>& m,
    const lduAddressing& addr,
    const lduInterfacePtrsList& interfaces
) const
{
    os  << "** Matrix **" << endl;
    const labelUList& upperAddr = addr.upperAddr();
    const labelUList& lowerAddr = addr.lowerAddr();
    const labelUList& ownerStartAddr = addr.ownerStartAddr();
    const labelUList& losortAddr = addr.losortAddr();

    const scalarField& lower = m.lower();
    const scalarField& upper = m.upper();
    const Field<Type>& source = m.source();
    const scalarField& diag = m.diag();


    // Invert patch addressing
    labelListList cellToPatch(addr.size());
    labelListList cellToPatchFace(addr.size());
    {
        forAll(interfaces, patchi)
        {
            if (interfaces.set(patchi))
            {
                const labelUList& fc = interfaces[patchi].faceCells();

                forAll(fc, i)
                {
                    cellToPatch[fc[i]].append(patchi);
                    cellToPatchFace[fc[i]].append(i);
                }
            }
        }
    }

    forAll(source, celli)
    {
        os  << "cell:" << celli << " diag:" << diag[celli]
            << " source:" << source[celli] << endl;

        label startLabel = ownerStartAddr[celli];
        label endLabel = ownerStartAddr[celli + 1];

        for (label facei = startLabel; facei < endLabel; facei++)
        {
            os  << "    face:" << facei
                << " nbr:" << upperAddr[facei] << " coeff:" << upper[facei]
                << endl;
        }

        startLabel = addr.losortStartAddr()[celli];
        endLabel = addr.losortStartAddr()[celli + 1];

        for (label i = startLabel; i < endLabel; i++)
        {
            label facei = losortAddr[i];
            os  << "    face:" << facei
                << " nbr:" << lowerAddr[facei] << " coeff:" << lower[facei]
                << endl;
        }

        forAll(cellToPatch[celli], i)
        {
            label patchi = cellToPatch[celli][i];
            label patchFacei = cellToPatchFace[celli][i];

            os  << "    patch:" << patchi
                << " patchface:" << patchFacei
                << " intcoeff:" << m.internalCoeffs()[patchi][patchFacei]
                << " boucoeff:" << m.boundaryCoeffs()[patchi][patchFacei]
                << endl;
        }
    }
    forAll(m.internalCoeffs(), patchi)
    {
        if (m.internalCoeffs().set(patchi))
        {
            const labelUList& fc = addr.patchAddr(patchi);

            os  << "patch:" << patchi
                //<< " type:" << interfaces[patchi].type()
                << " size:" << fc.size() << endl;
            if
            (
                interfaces.set(patchi)
             && isA<processorLduInterface>(interfaces[patchi])
            )
            {
                const processorLduInterface& ppp =
                    refCast<const processorLduInterface>(interfaces[patchi]);
                os  << "(processor with my:" << ppp.myProcNo()
                    << " nbr:" << ppp.neighbProcNo()
                    << ")" << endl;
            }

            forAll(fc, i)
            {
                os  << "    cell:" << fc[i]
                    << " int:" << m.internalCoeffs()[patchi][i]
                    << " boun:" << m.boundaryCoeffs()[patchi][i]
                    << endl;
            }
        }
    }


    lduInterfaceFieldPtrsList interfaceFields =
        m.psi().boundaryField().scalarInterfaces();
    forAll(interfaceFields, inti)
    {
        if (interfaceFields.set(inti))
        {
            os  << "interface:" << inti
                << " if type:" << interfaceFields[inti].interface().type()
                << " fld type:" << interfaceFields[inti].type() << endl;
        }
    }

    os  << "** End of Matrix **" << endl;
}


template<class GeoField>
void Foam::dynamicOversetFvMesh::correctCoupledBoundaryConditions(GeoField& fld)
{
    typename GeoField::Boundary& bfld = fld.boundaryFieldRef();

    const label nReq = Pstream::nRequests();

    forAll(bfld, patchi)
    {
        if (bfld[patchi].coupled())
        {
            //Pout<< "initEval of " << bfld[patchi].patch().name() << endl;
            bfld[patchi].initEvaluate(Pstream::defaultCommsType);
        }
    }

    // Block for any outstanding requests
    if (Pstream::parRun())
    {
        Pstream::waitRequests(nReq);
    }

    forAll(bfld, patchi)
    {
        if (bfld[patchi].coupled())
        {
            //Pout<< "eval of " << bfld[patchi].patch().name() << endl;
            bfld[patchi].evaluate(Pstream::defaultCommsType);
        }
    }
}


template<class GeoField>
void Foam::dynamicOversetFvMesh::checkCoupledBC(const GeoField& fld)
{
    Pout<< "** starting checking of " << fld.name() << endl;

    GeoField fldCorr(fld.name()+"_correct", fld);
    correctCoupledBoundaryConditions(fldCorr);

    const typename GeoField::Boundary& bfld = fld.boundaryField();
    const typename GeoField::Boundary& bfldCorr = fldCorr.boundaryField();

    forAll(bfld, patchi)
    {
        const typename GeoField::Patch& pfld = bfld[patchi];
        const typename GeoField::Patch& pfldCorr = bfldCorr[patchi];

        Pout<< "Patch:" << pfld.patch().name() << endl;

        forAll(pfld, i)
        {
            if (pfld[i] != pfldCorr[i])
            {
                Pout<< "    " << i << "  orig:" << pfld[i]
                    << " corrected:" << pfldCorr[i] << endl;
            }
        }
    }
    Pout<< "** end of " << fld.name() << endl;
}


// ************************************************************************* //

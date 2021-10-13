/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "dynamicOversetBlastFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "cellCellStencilObject.H"
#include "zeroGradientFvPatchFields.H"
#include "lduPrimitiveProcessorInterface.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicOversetBlastFvMesh, 0);
    addToRunTimeSelectionTable(dynamicBlastFvMesh, dynamicOversetBlastFvMesh, dictionary);
//     addToRunTimeSelectionTable(dynamicBlastFvMesh, dynamicOversetBlastFvMesh, doInit);
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicOversetBlastFvMesh::updateAddressing() const
{
    const cellCellStencilObject& overlap = Stencil::New(*this);

    // The (processor-local part of the) stencil determines the local
    // faces to add to the matrix. tbd: parallel
    const labelListList& stencil = overlap.cellStencil();

    // Get the base addressing
    //const lduAddressing& baseAddr = dynamicMotionSolverFvMesh::lduAddr();
    const lduAddressing& baseAddr = dynamicBlastFvMesh::lduAddr();

    // Add to the base addressing
    labelList lowerAddr;
    labelList upperAddr;
    label nExtraFaces;

    const globalIndex globalNumbering(baseAddr.size());
    labelListList localFaceCells;
    labelListList remoteFaceCells;

    labelList globalCellIDs(overlap.cellInterpolationMap().constructSize());
    forAll(baseAddr, cellI)
    {
        globalCellIDs[cellI] = globalNumbering.toGlobal(cellI);
    }
    overlap.cellInterpolationMap().distribute(globalCellIDs);


    reverseFaceMap_ = fvMeshPrimitiveLduAddressing::addAddressing
    (
        baseAddr,
        stencil,
        nExtraFaces,
        lowerAddr,
        upperAddr,
        stencilFaces_,
        globalNumbering,
        globalCellIDs,
        localFaceCells,
        remoteFaceCells
    );

    if (debug)
    {
        Pout<< "dynamicOversetBlastFvMesh::update() : extended addressing from"
            << " nFaces:" << baseAddr.lowerAddr().size()
            << " to nFaces:" << lowerAddr.size()
            << " nExtraFaces:" << nExtraFaces << endl;
    }


    // Now for the tricky bits. We want to hand out processor faces according
    // to the localFaceCells/remoteFaceCells. Ultimately we need
    // per entry in stencil:
    // - the patch (or -1 for internal faces)
    // - the face (is either an internal face index or a patch face index)

    stencilPatches_.setSize(stencilFaces_.size());

    // Per processor to owner (local)/neighbour (remote)
    List<DynamicList<label>> procOwner(Pstream::nProcs());
    List<DynamicList<label>> dynProcNeighbour(Pstream::nProcs());
    forAll(stencil, celli)
    {
        const labelList& nbrs = stencil[celli];
        stencilPatches_[celli].setSize(nbrs.size());
        stencilPatches_[celli] = -1;

        forAll(nbrs, nbri)
        {
            if (stencilFaces_[celli][nbri] == -1)
            {
                const label nbrCelli = nbrs[nbri];
                label globalNbr = globalCellIDs[nbrCelli];
                label proci = globalNumbering.whichProcID(globalNbr);
                label remoteCelli = globalNumbering.toLocal(proci, globalNbr);

                // Overwrite the face to be a patch face
                stencilFaces_[celli][nbri] = procOwner[proci].size();
                stencilPatches_[celli][nbri] = proci;
                procOwner[proci].append(celli);
                dynProcNeighbour[proci].append(remoteCelli);

                //Pout<< "From neighbour proc:" << proci
                //    << " allocating patchFace:" << stencilFaces_[celli][nbri]
                //    << " to get remote cell " << remoteCelli
                //    << " onto local cell " << celli << endl;
            }
        }
    }
    labelListList procNeighbour(dynProcNeighbour.size());
    forAll(procNeighbour, i)
    {
        procNeighbour[i] = std::move(dynProcNeighbour[i]);
    }
    labelListList mySendCells;
    Pstream::exchange<labelList, label>(procNeighbour, mySendCells);

    label nbri = 0;
    forAll(procOwner, proci)
    {
        if (procOwner[proci].size())
        {
            nbri++;
        }
        if (mySendCells[proci].size())
        {
            nbri++;
        }
    }
    remoteStencilInterfaces_.setSize(nbri);
    nbri = 0;

    // E.g. if proc1 needs some data from proc2 and proc2 needs some data from
    //      proc1. We first want the pair : proc1 receive and proc2 send
    //      and then the pair : proc1 send, proc2 receive


    labelList procToInterface(Pstream::nProcs(), -1);

    forAll(procOwner, proci)
    {
        if (proci < Pstream::myProcNo() && procOwner[proci].size())
        {
            if (debug)
            {
                Pout<< "Adding interface " << nbri
                    << " to receive my " << procOwner[proci].size()
                    << " from " << proci << endl;
            }
            procToInterface[proci] = nbri;
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    procOwner[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+2
                )
            );
        }
        else if (proci > Pstream::myProcNo() && mySendCells[proci].size())
        {
            if (debug)
            {
                Pout<< "Adding interface " << nbri
                    << " to send my " << mySendCells[proci].size()
                    << " to " << proci << endl;
            }
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    mySendCells[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+2
                )
            );
        }
    }
    forAll(procOwner, proci)
    {
        if (proci > Pstream::myProcNo() && procOwner[proci].size())
        {
            if (debug)
            {
                Pout<< "Adding interface " << nbri
                    << " to receive my " << procOwner[proci].size()
                    << " from " << proci << endl;
            }
            procToInterface[proci] = nbri;
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    procOwner[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+3
                )
            );
        }
        else if (proci < Pstream::myProcNo() && mySendCells[proci].size())
        {
            if (debug)
            {
                Pout<< "Adding interface " << nbri
                    << " to send my " << mySendCells[proci].size()
                    << " to " << proci << endl;
            }
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    mySendCells[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+3
                )
            );
        }
    }


    // Rewrite stencilPatches now we know the actual interface (procToInterface)
    for (auto& patches : stencilPatches_)
    {
        for (auto& interface : patches)
        {
            if (interface != -1)
            {
                interface = procToInterface[interface]+boundary().size();
            }
        }
    }


    // Get addressing and interfaces of all interfaces


    List<const labelUList*> patchAddr;
    {
        const fvBoundaryMesh& fvp = boundary();

        patchAddr.setSize(fvp.size() + remoteStencilInterfaces_.size());

        //allInterfaces_ = dynamicMotionSolverFvMesh::interfaces();
        allInterfaces_ = dynamicBlastFvMesh::interfaces();
        allInterfaces_.setSize(patchAddr.size());

        forAll(fvp, patchI)
        {
            patchAddr[patchI] = &fvp[patchI].faceCells();
        }
        forAll(remoteStencilInterfaces_, i)
        {
            label patchI = fvp.size()+i;
            const lduPrimitiveProcessorInterface& pp =
                remoteStencilInterfaces_[i];

            //Pout<< "at patch:" << patchI
            //    << " have procPatch:" << pp.type()
            //    << " from:" << pp.myProcNo()
            //    << " to:" << pp.neighbProcNo()
            //    << " with fc:" << pp.faceCells().size() << endl;

            patchAddr[patchI] = &pp.faceCells();
            allInterfaces_.set(patchI, &pp);
        }
    }
    const lduSchedule ps
    (
        lduPrimitiveMesh::nonBlockingSchedule<processorLduInterface>
        (
            allInterfaces_
        )
    );

    lduPtr_.reset
    (
        new fvMeshPrimitiveLduAddressing
        (
            nCells(),
            std::move(lowerAddr),
            std::move(upperAddr),
            patchAddr,
            ps
        )
   );


    // Check
    if (debug)
    {
        const lduAddressing& addr = lduPtr_();  //this->lduAddr();

        Pout<< "Adapted addressing:"
            << " lower:" << addr.lowerAddr().size()
            << " upper:" << addr.upperAddr().size() << endl;

        // Using lduAddressing::patch
        forAll(patchAddr, patchI)
        {
            Pout<< "    " << patchI << "\tpatchAddr:"
                << addr.patchAddr(patchI).size()
                << endl;
        }

        // Using interfaces
        const lduInterfacePtrsList& iFaces = allInterfaces_;
        Pout<< "Adapted interFaces:" << iFaces.size() << endl;
        forAll(iFaces, patchI)
        {
            if (iFaces.set(patchI))
            {
                Pout<< "    " << patchI << "\tinterface:"
                    << iFaces[patchI].type() << endl;
            }
        }
    }

    return true;
}


Foam::scalar Foam::dynamicOversetBlastFvMesh::cellAverage
(
    const labelList& types,
    const labelList& nbrTypes,
    const scalarField& norm,
    const scalarField& nbrNorm,
    const label celli,
    PackedBoolList& isFront
) const
{
    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();
    const cell& cFaces = cells()[celli];

    scalar avg = 0.0;
    label n = 0;
    label nFront = 0;
    for (const label facei : cFaces)
    {
        if (isInternalFace(facei))
        {
            label nbrCelli = (own[facei] == celli ? nei[facei] : own[facei]);
            if (norm[nbrCelli] == -GREAT)
            {
                // Invalid neighbour. Add to front
                if (isFront.set(facei))
                {
                    nFront++;
                }
            }
            else
            {
                // Valid neighbour. Add to average
                avg += norm[nbrCelli];
                n++;
            }
        }
        else
        {
            if (nbrNorm[facei-nInternalFaces()] == -GREAT)
            {
                if (isFront.set(facei))
                {
                    nFront++;
                }
            }
            else
            {
                avg += nbrNorm[facei-nInternalFaces()];
                n++;
            }
        }
    }

    if (n > 0)
    {
        return avg/n;
    }
    else
    {
        return norm[celli];
    }
}


void Foam::dynamicOversetBlastFvMesh::writeAgglomeration
(
    const GAMGAgglomeration& agglom
) const
{
    labelList cellToCoarse(identity(nCells()));
    labelListList coarseToCell(invertOneToMany(nCells(), cellToCoarse));

    // Write initial agglomeration
    {
        volScalarField scalarAgglomeration
        (
            IOobject
            (
                "agglomeration",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimensionedScalar(dimless, Zero)
        );
        scalarField& fld = scalarAgglomeration.primitiveFieldRef();
        forAll(fld, celli)
        {
            fld[celli] = cellToCoarse[celli];
        }
        fld /= max(fld);
        correctBoundaryConditions
        <
            volScalarField,
            oversetFvPatchField<scalar>
        >(scalarAgglomeration.boundaryFieldRef(), false);
        scalarAgglomeration.write();

        Info<< "Writing initial cell distribution to "
            << this->time().timeName() << endl;
    }


    for (label level = 0; level < agglom.size(); level++)
    {
        const labelList& addr = agglom.restrictAddressing(level);
        label coarseSize = max(addr)+1;

        Info<< "Level : " << level << endl
            << returnReduce(addr.size(), sumOp<label>()) << endl
            << "    current size      : "
            << returnReduce(addr.size(), sumOp<label>()) << endl
            << "    agglomerated size : "
            << returnReduce(coarseSize, sumOp<label>()) << endl;

        forAll(addr, fineI)
        {
            const labelList& cellLabels = coarseToCell[fineI];
            forAll(cellLabels, i)
            {
                cellToCoarse[cellLabels[i]] = addr[fineI];
            }
        }
        coarseToCell = invertOneToMany(coarseSize, cellToCoarse);

        // Write agglomeration
        {
            volScalarField scalarAgglomeration
            (
                IOobject
                (
                    "agglomeration_" + Foam::name(level),
                    this->time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                *this,
                dimensionedScalar(dimless, Zero)
            );
            scalarField& fld = scalarAgglomeration.primitiveFieldRef();
            forAll(fld, celli)
            {
                fld[celli] = cellToCoarse[celli];
            }
            //if (normalise)
            //{
            //    fld /= max(fld);
            //}
            correctBoundaryConditions
            <
                volScalarField,
                oversetFvPatchField<scalar>
            >(scalarAgglomeration.boundaryFieldRef(), false);
            scalarAgglomeration.write();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicOversetBlastFvMesh::dynamicOversetBlastFvMesh
(
    const IOobject& io//,
//     const bool doInit
)
:
    dynamicMotionSolverBlastFvMesh(io)//, doInit)
{
//     if (doInit)
    {
//         init(false);    // do not initialise lower levels
        init(true);
    }
}


bool Foam::dynamicOversetBlastFvMesh::init(const bool doInit)
{
//     if (doInit)
//     {
//         dynamicMotionSolverBlastFvMesh::init(doInit);
//     }

    active_ = false;

    // Load stencil (but do not update)
    (void) Stencil::New(dynamicCast<const fvMesh&>(*this), false);

    // Assume something changed
    return true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicOversetBlastFvMesh::~dynamicOversetBlastFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::lduAddressing& Foam::dynamicOversetBlastFvMesh::lduAddr() const
{
    if (!active_)
    {
        //return dynamicMotionSolverBlastFvMesh::lduAddr();
        return dynamicBlastFvMesh::lduAddr();
    }
    if (!lduPtr_.valid())
    {
        // Build extended addressing
        updateAddressing();
    }
    return lduPtr_();
}


Foam::lduInterfacePtrsList Foam::dynamicOversetBlastFvMesh::interfaces() const
{
    if (!active_)
    {
        //return dynamicMotionSolverBlastFvMesh::interfaces();
        return dynamicBlastFvMesh::interfaces();
    }
    if (!lduPtr_.valid())
    {
        // Build extended addressing
        updateAddressing();
    }
    return allInterfaces_;
}


const Foam::fvMeshPrimitiveLduAddressing&
Foam::dynamicOversetBlastFvMesh::primitiveLduAddr() const
{
    if (!lduPtr_.valid())
    {
        FatalErrorInFunction
            << "Extended addressing not allocated" << abort(FatalError);
    }

    return lduPtr_();
}


bool Foam::dynamicOversetBlastFvMesh::update()
{
    //if (dynamicMotionSolverBlastFvMesh::update())
    if (dynamicMotionSolverBlastFvMesh::update())
    {
        // Calculate the local extra faces for the interpolation. Note: could
        // let demand-driven lduAddr() trigger it but just to make sure.
        updateAddressing();

        // Addressing and/or weights have changed. Make interpolated cells
        // up to date with donors
        interpolateFields();

        return true;
    }

    return false;
}


Foam::word Foam::dynamicOversetBlastFvMesh::baseName(const word& name)
{
    if (label(name.find("_0", name.size() - 2)) >= 0)
    {
        return baseName(name.substr(0, name.size()-2));
    }

    return name;
}


bool Foam::dynamicOversetBlastFvMesh::interpolateFields()
{
    // Add the stencil suppression list
    wordHashSet suppressed(Stencil::New(*this).nonInterpolatedFields());

    // Use whatever the solver has set up as suppression list
    const dictionary* dictPtr
    (
        this->schemesDict().subDictPtr("oversetInterpolationSuppressed")
    );
    if (dictPtr)
    {
        suppressed.insert(dictPtr->toc());
    }

    interpolate<volScalarField>(suppressed);
    interpolate<volVectorField>(suppressed);
    interpolate<volSphericalTensorField>(suppressed);
    interpolate<volSymmTensorField>(suppressed);
    interpolate<volTensorField>(suppressed);

    return true;
}



bool Foam::dynamicOversetBlastFvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    bool ok = dynamicMotionSolverBlastFvMesh::writeObject(fmt, ver, cmp, write);

    // For postprocessing : write cellTypes and zoneID
    {
        const cellCellStencilObject& overlap = Stencil::New(*this);

        const labelUList& cellTypes = overlap.cellTypes();

        volScalarField volTypes
        (
            IOobject
            (
                "cellTypes",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimensionedScalar(dimless, Zero),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(volTypes.internalField(), cellI)
        {
            volTypes[cellI] = cellTypes[cellI];
        }
        volTypes.correctBoundaryConditions();
        volTypes.writeObject(fmt, ver, cmp, write);
    }
    {
        volScalarField volZoneID
        (
            IOobject
            (
                "zoneID",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimensionedScalar(dimless, Zero),
            zeroGradientFvPatchScalarField::typeName
        );

        const cellCellStencilObject& overlap = Stencil::New(*this);
        const labelIOList& zoneID = overlap.zoneID();

        forAll(zoneID, cellI)
        {
            volZoneID[cellI] = zoneID[cellI];
        }
        volZoneID.correctBoundaryConditions();
        volZoneID.writeObject(fmt, ver, cmp, write);
    }
    if (debug)
    {
        const cellCellStencilObject& overlap = Stencil::New(*this);
        const labelIOList& zoneID = overlap.zoneID();
        const labelListList& cellStencil = overlap.cellStencil();

        // Get remote zones
        labelList donorZoneID(zoneID);
        overlap.cellInterpolationMap().distribute(donorZoneID);

        // Get remote cellCentres
        pointField cc(C());
        overlap.cellInterpolationMap().distribute(cc);

        volScalarField volDonorZoneID
        (
            IOobject
            (
                "donorZoneID",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimensionedScalar("minOne", dimless, scalar(-1)),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(cellStencil, cellI)
        {
            const labelList& stencil = cellStencil[cellI];
            if (stencil.size())
            {
                volDonorZoneID[cellI] = donorZoneID[stencil[0]];
                for (label i = 1; i < stencil.size(); i++)
                {
                    if (donorZoneID[stencil[i]] != volDonorZoneID[cellI])
                    {
                        WarningInFunction << "Mixed donor meshes for cell "
                            << cellI << " at " << C()[cellI]
                            << " donors:" << UIndirectList<point>(cc, stencil)
                            << endl;
                        volDonorZoneID[cellI] = -2;
                    }
                }
            }
        }
        //- Do not correctBoundaryConditions since re-interpolates!
        //volDonorZoneID.correctBoundaryConditions();
        volDonorZoneID.writeObject(fmt, ver, cmp, write);
    }

    return ok;
}


// ************************************************************************* //

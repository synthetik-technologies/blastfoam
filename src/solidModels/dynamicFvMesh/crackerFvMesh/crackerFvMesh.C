/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "crackerFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "mapPolyMesh.H"
//#include "materialInterface.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(crackerFvMesh, 0);
    addToRunTimeSelectionTable(topoChangerFvMesh, crackerFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::crackerFvMesh::removeZones()
{
    if (debug)
    {
        Info<< "void polyMesh::removeZones(): "
            << "Removing point, face and cell zones."
            << endl;
    }

    // Remove the zones and reset writing
    meshPointZones& pZones =
        const_cast<meshPointZones&>(this->pointZones());
    pZones.clear();
    pZones.setSize(0);
    pZones.writeOpt() = IOobject::NO_WRITE;

    meshFaceZones& fZones =
        const_cast<meshFaceZones&>(this->faceZones());
    fZones.clear();
    fZones.setSize(0);
    fZones.writeOpt() = IOobject::NO_WRITE;

    meshCellZones& cZones =
        const_cast<meshCellZones&>(this->cellZones());
    cZones.clear();
    cZones.setSize(0);
    cZones.writeOpt() = IOobject::NO_WRITE;

    polyMesh::clearOut();
}


void Foam::crackerFvMesh::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

    Info<< "Time = " << time().timeName() << endl
        << "Adding topo modifier to the mesh" << endl;

    const word crackPatchName(dict_.lookup("crackPatch"));
    const label crackPatchIndex = boundaryMesh().findPatchID(crackPatchName);

    if (crackPatchIndex < 0)
    {
        FatalErrorInFunction
            << "Crack patch not found in boundary"
            << abort(FatalError);
    }

    // Add zones
    if (faceZones().findZoneID(crackPatchName + "Zone") == -1)
    {
        Info << "Adding the crack faceZone" << endl;

        // Copy points zones from the mesh

        List<pointZone*> pz(pointZones().size());

        forAll(pz, zoneI)
        {
            pz[zoneI] =
                new pointZone
                (
                    pointZones()[zoneI],
                    pointZones()[zoneI],
                    zoneI,
                    pointZones()
                );
        }

        // Copy face zones from the mesh and add a crack zone at the end

        List<faceZone*> fz(faceZones().size() + 1);

        for (label zoneI = 0; zoneI < fz.size() - 1; zoneI++)
        {
            fz[zoneI] =
                new faceZone
                (
                    faceZones()[zoneI],
                    faceZones()[zoneI],
                    faceZones()[zoneI].flipMap(),
                    zoneI,
                    faceZones()
                );
        }

        // Add crack face zone at the end
        fz[fz.size() - 1] =
            new faceZone
            (
                crackPatchName + "Zone",
                labelList(0),
                boolList(0),
                fz.size() - 1,
                faceZones()
            );

        // Copy cell zones from the mesh

        List<cellZone*> cz(cellZones().size());

        forAll(cz, zoneI)
        {
            cz[zoneI] =
                new cellZone
                (
                    cellZones()[zoneI],
                    cellZones()[zoneI],
                    zoneI,
                    cellZones()
                );
        }

        // Remove previous zones
        removeZones();

        // Add the zones to the mesh
        addZones(pz, fz, cz);
    }
    else
    {
        Info << "Face zones already present" << endl;
    }

    // Add a topology modifier
    if (topoChanger_.size() == 0)
    {
        Info << "Adding topology modifiers" << endl;
        topoChanger_.setSize(1);
        topoChanger_.set
        (
            0,
            new faceCracker
            (
                "cracker",
                0,
                topoChanger_,
                crackPatchName + "Zone",
                crackPatchName //,
                //openPatchName
            )
        );

        topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    }
    else
    {
        Info<< "void crackerFvMesh::addZonesAndModifiers() : "
            << "Modifiers already present."
            << endl;
    }

    // Write mesh
    write();
}

void Foam::crackerFvMesh::makeRegions() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (regionsPtr_)
    {
        FatalErrorInFunction
            << "regions already exist"
            << abort(FatalError);
    }

    regionsPtr_ = new regionSplit(*this);
}

void Foam::crackerFvMesh::makeNCellsInRegion() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (nCellsInRegionPtr_)
    {
        FatalErrorInFunction
            << "number of cells in regions already exist"
            << abort(FatalError);
    }

    nCellsInRegionPtr_ = new labelList(regions().nRegions(), 0);

    labelList& nCellsInRegion = *nCellsInRegionPtr_;

    const labelList& regs = regions();

    forAll(regs, cellI)
    {
        nCellsInRegion[regs[cellI]]++;
    }
}

void Foam::crackerFvMesh::makeGlobalCrackFaceCentresAndSizes() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (globalCrackFaceCentresPtr_ || globalCrackFaceSizesPtr_)
    {
        FatalErrorInFunction
            << "global crack face centres and sizes already exist"
            << abort(FatalError);
    }


    // Number of faces in global crack
    labelList sizes(Pstream::nProcs(), 0);
    sizes[Pstream::myProcNo()] = boundaryMesh()[crackPatchID_.index()].size();

    if (Pstream::parRun())
    {
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::commsTypes::blocking,
                        procI,
                        sizeof(label)
                    );

                    toProc << sizes[Pstream::myProcNo()];
                }
            }
        }

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::commsTypes::blocking,
                        procI,
                        sizeof(label)
                    );

                    fromProc >> sizes[procI];
                }
            }
        }
    }

    label globalCrackSize = sum(sizes);

    globalCrackFaceCentresPtr_ =
        new vectorField(globalCrackSize, vector::zero);
    vectorField& globalCrackFaceCentres = *globalCrackFaceCentresPtr_;

    globalCrackFaceSizesPtr_ =
        new scalarField(globalCrackSize, 0);
    scalarField& globalCrackFaceSizes = *globalCrackFaceSizesPtr_;

    localCrackStart_ = 0;
    for (label procI = 0; procI < Pstream::myProcNo(); procI++)
    {
        localCrackStart_ += sizes[procI];
    }

//     const vectorField& crackCf =
//         boundaryMesh()[crackPatchID_.index()].faceCentres();
    const vectorField crackCf
    (
        boundaryMesh()[crackPatchID_.index()].faceCentres()
    );

    // Calc face sizes
//     const vectorField& crackSf =
//         boundaryMesh()[crackPatchID_.index()].faceAreas();
    const vectorField crackSf
    (
        boundaryMesh()[crackPatchID_.index()].faceAreas()
    );

    scalarField delta(crackSf.size(), 0);

    if (nGeometricD() == 3)
    {
        delta = Foam::sqrt(mag(crackSf));
    }
    else
    {
        scalar thickness = 0.0;
        const Vector<label>& directions = geometricD();
        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = bounds().span()[dir];
                break;
            }
        }

        delta = mag(crackSf)/thickness;
    }

    label j=0;
    for
    (
        label i=localCrackStart_;
        i<(localCrackStart_ + sizes[Pstream::myProcNo()]);
        i++
    )
    {
        globalCrackFaceCentres[i] = crackCf[j];
        globalCrackFaceSizes[i] = delta[j];
        j++;
    }

    // Parallel data exchange: collect crack face centres and sizes
    // on all processors
    reduce(globalCrackFaceCentres, sumOp<List<vector>>());
    reduce(globalCrackFaceSizes, sumOp<List<scalar>>());
}


void Foam::crackerFvMesh::makeGlobalCrackFaceAddressing() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (globalCrackFaceAddressingPtr_)
    {
        FatalErrorInFunction
            << "global crack face addressing already exists"
            << abort(FatalError);
    }

    const vectorField& gcfc = globalCrackFaceCentres();
    const scalarField& gcfs = globalCrackFaceSizes();

    globalCrackFaceAddressingPtr_ = new labelList(gcfc.size(), -1);
    labelList& gcfa = *globalCrackFaceAddressingPtr_;

    forAll(gcfa, faceI)
    {
        if (gcfa[faceI] < 0)
        {
            forAll(gcfc, fI)
            {
                if ((fI != faceI) && (gcfa[fI] < 0))
                {
                    if (mag(gcfc[faceI] - gcfc[fI]) < 1e-3*gcfs[faceI])
                    {
                        gcfa[faceI] = fI;
                        gcfa[fI] = faceI;
                        break;
                    }
                }
            }
        }
    }

    // Check addressing
    forAll(gcfa, faceI)
    {
        if (gcfa[faceI] < 0)
        {
            FatalErrorInFunction
            << "problem with defining global crack face addressing"
            << abort(FatalError);
        }
    }
}


void Foam::crackerFvMesh::makeFaceBreakerLaw() const
{
    if (lawPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    lawPtr_ =
        faceBreakerLaw::New
        (
            "law",
            *this,
            dict_.subDict("law")
        );
}


void Foam::crackerFvMesh::perturbFieldOnNewCrackFaces
(
    const labelList& faceMap,
    const labelList& facesToBreak,
    const labelList& coupledFacesToBreak,
    const word& fieldName
) const
{
    // Create a hash set from the list of faces for easy searching
    const labelHashSet facesToBreakSet(facesToBreak);
    const labelHashSet coupledFacesToBreakSet(coupledFacesToBreak);

    // Cast the mesh to a crackerFvMesh

    const crackerFvMesh& mesh = refCast<const crackerFvMesh>(*this);

    // Lookup field from object registry
    volVectorField& field
    (
        const_cast<volVectorField&>
        (
            mesh.thisDb().objectRegistry::lookupObject<volVectorField>
            (
                fieldName
            )
        )
    );

    const label cohesivePatchID = crackPatchID_.index();
    const label start = mesh.boundaryMesh()[cohesivePatchID].start();
    const label cohesivePatchSize = mesh.boundaryMesh()[cohesivePatchID].size();

    Info<< "    Perturbing " << fieldName << " on new crack faces" << endl;

    // Local crack field
    Field<vector> fieldpI
    (
        field.boundaryField()[cohesivePatchID].patchInternalField()
    );

    // Global crack fields
    Field<vector> gFieldpI(mesh.globalCrackField(fieldpI));

    //const labelList& gcfa = mesh.globalCrackFaceAddressing();

    label globalIndex = mesh.localCrackStart();

    for (label i = 0; i < cohesivePatchSize; i++)
    {
        label oldFaceIndex = faceMap[start + i];

        if
        (
            facesToBreakSet.found(oldFaceIndex)
         || coupledFacesToBreakSet.found(oldFaceIndex)
        )
        {
            // We would like the crack to be slightly open so as any fluid
            // mesh cells grown inside do not have zero volume when created
            // So we will slightly perturb the displacement field on the
            // new crack face
            // We will add a small displacement in the negative face normal
            // direction
            const vector& faceN =
                mesh.boundaryMesh()[cohesivePatchID].faceNormals()[i];

            field.boundaryFieldRef()[cohesivePatchID][i] -= 1e-12*faceN;

            globalIndex++;
        }
        else
        {
            globalIndex++;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::crackerFvMesh::crackerFvMesh
(
    const IOobject& io
)
:
    topoChangerFvMesh(io),
    dict_(dynamicMeshDict().optionalSubDict(type() + "Coeffs")),
    topoChangeMap_(),
    crackPatchID_
    (
        dict_.lookup<word>("crackPatch"),
        boundaryMesh()
    ),
    lawPtr_(NULL),
    regionsPtr_(NULL),
    nCellsInRegionPtr_(NULL),
    globalCrackFaceCentresPtr_(NULL),
    globalCrackFaceSizesPtr_(NULL),
    localCrackStart_(-1),
    globalCrackFaceAddressingPtr_(NULL)
{
    // Add zones and mesh modifiers
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::crackerFvMesh::~crackerFvMesh()
{
    deleteDemandDrivenData(regionsPtr_);
    deleteDemandDrivenData(nCellsInRegionPtr_);
    deleteDemandDrivenData(globalCrackFaceCentresPtr_);
    deleteDemandDrivenData(globalCrackFaceSizesPtr_);
    deleteDemandDrivenData(globalCrackFaceAddressingPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::crackerFvMesh::setBreak
(
    const labelList& facesToBreak,
    const boolList& faceFlip,
    const labelList& coupledFacesToBreak
)
{
    faceCracker& fc = refCast<faceCracker>(topoChanger_[0]);

    fc.setBreak(facesToBreak, faceFlip, coupledFacesToBreak);
}


bool Foam::crackerFvMesh::update()
{
    // Clearout the law demand driven data
    faceBreaker().clearOut();

    // Get faces to break from the law
    const labelList& facesToBreak = faceBreaker().facesToBreak();
    const boolList facesToBreakFlip = boolList(facesToBreak.size(), false);
    const labelList& coupledFacesToBreak = faceBreaker().coupledFacesToBreak();

    // All processors must know if a topological change will occur
    label nFacesToBreak = facesToBreak.size();
    label nCoupledFacesToBreak = coupledFacesToBreak.size();
    reduce(nFacesToBreak, maxOp<label>());
    reduce(nCoupledFacesToBreak, maxOp<label>());

    if (nFacesToBreak || nCoupledFacesToBreak)
    {
        if (debug)
        {
            Pout<< "nFacesToBreak: " << nFacesToBreak << nl
                << "nCoupledFacesToBreak: " << nCoupledFacesToBreak << nl
                << "facesToBreak: " << facesToBreak << nl
                << "coupledFacesToBreak: " << coupledFacesToBreak << endl;
        }

        // Set faces to break
        setBreak(facesToBreak, facesToBreakFlip, coupledFacesToBreak);

        // Perform mesh topological change to break the faces
        topoChangeMap_ = topoChanger_.changeMesh(true);

        if (topoChangeMap_.valid())
        {
            deleteDemandDrivenData(regionsPtr_);
            deleteDemandDrivenData(nCellsInRegionPtr_);
            deleteDemandDrivenData(globalCrackFaceCentresPtr_);
            deleteDemandDrivenData(globalCrackFaceSizesPtr_);
            localCrackStart_ = -1;
            deleteDemandDrivenData(globalCrackFaceAddressingPtr_);
        }

        // Update field values on the new crack faces

        const labelList& faceMap = topoChangeMap().faceMap();

        Info<< "Updating field values on newly broken faces" << endl;

        updateVolFieldBrokenFaces<scalar>
        (
            faceMap, facesToBreak, coupledFacesToBreak
        );
        updateVolFieldBrokenFaces<vector>
        (
            faceMap, facesToBreak, coupledFacesToBreak
        );
        updateVolFieldBrokenFaces<tensor>
        (
            faceMap, facesToBreak, coupledFacesToBreak
        );
        updateVolFieldBrokenFaces<symmTensor>
        (
            faceMap, facesToBreak, coupledFacesToBreak
        );
//         updateVolFieldBrokenFaces<diagTensor>
//         (
//             faceMap, facesToBreak, coupledFacesToBreak
//         );
        updateVolFieldBrokenFaces<sphericalTensor>
        (
            faceMap, facesToBreak, coupledFacesToBreak
        );

        // Point fields should be recalculated in the solver

        // Clearout mechanical properties within interface corrector so they
        // will be regenerated
        // if (foundObject<materialInterface>("materialInterface"))
        // {
        //     materialInterface& interface =
        //         const_cast<materialInterface&>
        //         (
        //             lookupObject<materialInterface>("materialInterface")
        //         );

        //     interface.clearOut();
        // }

        // Note: after cracking, dead cell regions (small groups of cells
        // unconnected to the main mesh may have been created); it may be
        // required to set a reference to these dead cells in the solver after
        // or alternatively to delete them.
    }
    return bool(nFacesToBreak || nCoupledFacesToBreak);
}


const Foam::regionSplit& Foam::crackerFvMesh::regions() const
{
    if (!regionsPtr_)
    {
        makeRegions();
    }

    return *regionsPtr_;
}


Foam::label Foam::crackerFvMesh::nCellsInRegion(label regI) const
{
    if (!nCellsInRegionPtr_)
    {
        makeNCellsInRegion();
    }

    if ((regI < 0) || (regI >= regions().nRegions()))
    {
        FatalErrorInFunction
            << "region index is out of range"
            << abort(FatalError);
    }

    return (*nCellsInRegionPtr_)[regI];
}


const Foam::mapPolyMesh& Foam::crackerFvMesh::topoChangeMap() const
{
    return topoChangeMap_();
}


const Foam::vectorField& Foam::crackerFvMesh::globalCrackFaceCentres() const
{
    if (!globalCrackFaceCentresPtr_)
    {
        makeGlobalCrackFaceCentresAndSizes();
    }

    return *globalCrackFaceCentresPtr_;
}


const Foam::scalarField& Foam::crackerFvMesh::globalCrackFaceSizes() const
{
    if (!globalCrackFaceSizesPtr_)
    {
        makeGlobalCrackFaceCentresAndSizes();
    }

    return *globalCrackFaceSizesPtr_;
}


const Foam::labelList& Foam::crackerFvMesh::globalCrackFaceAddressing() const
{
    if (!globalCrackFaceAddressingPtr_)
    {
        makeGlobalCrackFaceAddressing();
    }

    return *globalCrackFaceAddressingPtr_;
}


Foam::label Foam::crackerFvMesh::localCrackStart() const
{
    if (localCrackStart_ == -1)
    {
        makeGlobalCrackFaceCentresAndSizes();
    }

    return localCrackStart_;
}


Foam::label Foam::crackerFvMesh::globalCrackSize() const
{
    return globalCrackFaceCentres().size();
}


const Foam::faceBreakerLaw& Foam::crackerFvMesh::faceBreaker() const
{
    if (!lawPtr_.valid())
    {
        makeFaceBreakerLaw();
    }

    return lawPtr_();
}


Foam::faceBreakerLaw& Foam::crackerFvMesh::faceBreaker()
{
    if (!lawPtr_.valid())
    {
        makeFaceBreakerLaw();
    }

    return lawPtr_();
}


// ************************************************************************* //

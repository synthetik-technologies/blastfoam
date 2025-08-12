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

#include "crackerFvMeshTopoChanger.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(cracker, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, cracker, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshTopoChangers::cracker::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action
    DebugInfo<< "Adding topo modifier to the mesh" << endl;

    const label crackPatchIndex =
        mesh().boundaryMesh().findIndex(crackPatch_);

    if (crackPatchIndex < 0)
    {
        FatalErrorInFunction
            << "Crack patch not found in boundary"
            << abort(FatalError);
    }

    // Add zones
    if (mesh().faceZones().findIndex(crackPatch_ + "Zone") == -1)
    {
        DebugInfo<< "Adding the crack faceZone" << endl;

        // Copy points zones from the mesh

        const pointZoneList& pzs = mesh().pointZones();
        List<pointZone*> newPzs(pzs.size());
        forAll(pzs, zoneI)
        {
            newPzs[zoneI] = pzs[zoneI].clone(pzs).ptr();
        }

        // Copy face zones from the mesh and add a crack zone at the end

        const faceZoneList& fzs = mesh().faceZones();
        List<faceZone*> newFzs(fzs.size() + 1);
        forAll(fzs, zoneI)
        {
            newFzs[zoneI] = fzs[zoneI].clone(fzs).ptr();
        }

        // Add crack face zone at the end
        newFzs[newFzs.size() - 1] =
            new faceZone
            (
                crackPatch_ + "Zone",
                labelList(0),
                boolList(0),
                fzs
            );

        // Copy cell zones from the mesh

        const cellZoneList& czs = mesh().cellZones();
        List<cellZone*> newCzs(czs.size());
        forAll(czs, zoneI)
        {
            newCzs[zoneI] = czs[zoneI].clone(czs).ptr();
        }

        // Remove previous zones
        mesh().pointZones().clear();
        mesh().faceZones().clear();
        mesh().cellZones().clear();

        // Add the zones to the mesh
        mesh().addZones(newPzs, newFzs, newCzs);
    }
    else
    {
        DebugInfo<< "Face zones already present" << endl;
    }

    // Add a topology modifier
    if (!topoChanger_.valid())
    {
        DebugInfo<< "Adding topology modifiers" << endl;
        topoChanger_.set
        (
            new faceCracker
            (
                mesh(),
                word(crackPatch_ + "Zone"),
                crackPatch_
            )
        );
    }
    else
    {
        DebugInfo<< "Modifiers already present." << endl;
    }

    // Write mesh
    mesh().write();
}

void Foam::fvMeshTopoChangers::cracker::makeRegions() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (regionsPtr_)
    {
        FatalErrorInFunction
            << "regions already exist"
            << abort(FatalError);
    }

    regionsPtr_ = new regionSplit(mesh());
}

void Foam::fvMeshTopoChangers::cracker::makeNCellsInRegion() const
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

void
Foam::fvMeshTopoChangers::cracker::makeGlobalCrackFaceCentresAndSizes() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (globalCrackFaceCentresPtr_ || globalCrackFaceSizesPtr_)
    {
        FatalErrorInFunction
            << "global crack face centres and sizes already exist"
            << abort(FatalError);
    }


    // Crack patch
    const label crackPatchID = mesh().boundaryMesh().findIndex(crackPatch_);
    const polyPatch& crackPatch = mesh().boundaryMesh()[crackPatchID];

    // Number of faces in global crack
    labelList sizes(Pstream::nProcs(), 0);
    sizes[Pstream::myProcNo()] = crackPatch.size();
    Pstream::gatherList(sizes);
    Pstream::scatterList(sizes);

    const label globalCrackSize = sum(sizes);

    globalCrackFaceCentresPtr_ = new vectorField(globalCrackSize, Zero);
    vectorField& crackFaceCentres = *globalCrackFaceCentresPtr_;

    globalCrackFaceSizesPtr_ = new scalarField(globalCrackSize, 0.0);
    scalarField& crackFaceSizes = *globalCrackFaceSizesPtr_;

    localCrackStart_ = 0;
    for (label procI = 0; procI < Pstream::myProcNo(); procI++)
    {
        localCrackStart_ += sizes[procI];
    }

    const vectorField::subField crackCf =
        mesh().boundaryMesh()[crackPatchID].faceCentres();

    // Calc face sizes
    const vectorField::subField crackSf =
        mesh().boundaryMesh()[crackPatchID].faceAreas();

    scalarField delta(crackSf.size(), 0.0);
    if (mesh().nGeometricD() == 3)
    {
        delta = Foam::sqrt(mag(crackSf));
    }
    else
    {
        scalar thickness = 0.0;
        const Vector<label>& directions = mesh().geometricD();
        const vector span = mesh().bounds().span();
        for (direction dir = 0; dir < vector::nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = span[dir];
                break;
            }
        }

        delta = mag(crackSf)/thickness;
    }


    if (Pstream::parRun())
    {
        List<vectorField> globalCrackFaceCentres(Pstream::nProcs());
        globalCrackFaceCentres[Pstream::myProcNo()] = crackCf;

        List<scalarField> globalCrackFaceSizes(Pstream::nProcs());
        globalCrackFaceSizes[Pstream::myProcNo()] = delta;

        Pstream::gatherList(globalCrackFaceCentres);
        Pstream::scatterList(globalCrackFaceCentres);
        Pstream::gatherList(globalCrackFaceSizes);
        Pstream::scatterList(globalCrackFaceSizes);

        label facei = 0;
        forAll(globalCrackFaceCentres, proci)
        {
            const vectorField& gCrackFaceCentres =
                globalCrackFaceCentres[proci];
            const scalarField& gCrackFaceSizes =
                globalCrackFaceSizes[proci];
            forAll(gCrackFaceCentres, fi)
            {
                crackFaceCentres[facei] = gCrackFaceCentres[fi];
                crackFaceSizes[facei] = gCrackFaceSizes[fi];
                facei++;
            }
        }
    }
    else
    {
        crackFaceCentres = crackCf;
        crackFaceSizes = delta;
    }
}


void Foam::fvMeshTopoChangers::cracker::makeGlobalCrackFaceAddressing() const
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


void Foam::fvMeshTopoChangers::cracker::perturbFieldOnNewCrackFaces
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

    // Cast the mesh to a cracker

    const fvMesh& mesh = this->mesh();

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

    const polyPatch& crackPatch = mesh().boundaryMesh()[crackPatch_];
    const label crackPatchID = crackPatch.index();
    const label start = crackPatch.start();

    Info<< "    Perturbing " << fieldName << " on new crack faces" << endl;

    // Local crack field
    Field<vector> fieldpI
    (
        field.boundaryField()[crackPatchID].patchInternalField()
    );

    // Global crack fields
    Field<vector> gFieldpI(this->globalCrackField(fieldpI));

    volVectorField::Patch& pf =
        field.boundaryFieldRefNoStoreOldTimes()[crackPatchID];

    forAll(crackPatch, fi)
    {
        label oldFaceIndex = faceMap[start + fi];

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
            pf[fi] -= 1e-12*crackPatch.faceNormals()[fi];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::cracker::cracker
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    regIOobject
    (
        IOobject
        (
            type(),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    ),
    fvMeshTopoChanger(mesh),
    dict_(dict),
    crackPatch_(dict.lookup<word>("crackPatch")),
    lawPtr_
    (
        faceBreakerLaw::New
        (
            "law",
            mesh,
            dict.optionalSubDict(faceBreakerLaw::typeName + "Coeffs")
        )
    ),
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

Foam::fvMeshTopoChangers::cracker::~cracker()
{
    deleteDemandDrivenData(regionsPtr_);
    deleteDemandDrivenData(nCellsInRegionPtr_);
    deleteDemandDrivenData(globalCrackFaceCentresPtr_);
    deleteDemandDrivenData(globalCrackFaceSizesPtr_);
    deleteDemandDrivenData(globalCrackFaceAddressingPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshTopoChangers::cracker::setBreak
(
    const labelList& facesToBreak,
    const boolList& faceFlip,
    const labelList& coupledFacesToBreak
)
{

    topoChanger_->setBreak(mesh(), facesToBreak, faceFlip, coupledFacesToBreak);
}


bool Foam::fvMeshTopoChangers::cracker::update()
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
        polyTopoChange meshMod(mesh());
        topoChanger_->setRefinement(meshMod);

        autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh());
        mesh().topoChange(map);

        {
            deleteDemandDrivenData(regionsPtr_);
            deleteDemandDrivenData(nCellsInRegionPtr_);
            deleteDemandDrivenData(globalCrackFaceCentresPtr_);
            deleteDemandDrivenData(globalCrackFaceSizesPtr_);
            deleteDemandDrivenData(globalCrackFaceAddressingPtr_);
        }

        // Update field values on the new crack faces

        const labelList& faceMap = map->faceMap();

        DebugInfo<< "Updating field values on newly broken faces" << endl;

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
    return
        returnReduce
        (
            bool(nFacesToBreak || nCoupledFacesToBreak),
            orOp<bool>()
        );
}


const Foam::regionSplit& Foam::fvMeshTopoChangers::cracker::regions() const
{
    if (!regionsPtr_)
    {
        makeRegions();
    }

    return *regionsPtr_;
}


Foam::label Foam::fvMeshTopoChangers::cracker::nCellsInRegion(label regI) const
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


const Foam::vectorField&
Foam::fvMeshTopoChangers::cracker::globalCrackFaceCentres() const
{
    if (!globalCrackFaceCentresPtr_)
    {
        makeGlobalCrackFaceCentresAndSizes();
    }

    return *globalCrackFaceCentresPtr_;
}


const Foam::scalarField&
Foam::fvMeshTopoChangers::cracker::globalCrackFaceSizes() const
{
    if (!globalCrackFaceSizesPtr_)
    {
        makeGlobalCrackFaceCentresAndSizes();
    }

    return *globalCrackFaceSizesPtr_;
}


const Foam::labelList&
Foam::fvMeshTopoChangers::cracker::globalCrackFaceAddressing() const
{
    if (!globalCrackFaceAddressingPtr_)
    {
        makeGlobalCrackFaceAddressing();
    }

    return *globalCrackFaceAddressingPtr_;
}


Foam::label
Foam::fvMeshTopoChangers::cracker::localCrackStart() const
{
    if (localCrackStart_ == -1)
    {
        makeGlobalCrackFaceCentresAndSizes();
    }

    return localCrackStart_;
}


Foam::label
Foam::fvMeshTopoChangers::cracker::globalCrackSize() const
{
    return globalCrackFaceCentres().size();
}


const Foam::faceBreakerLaw&
Foam::fvMeshTopoChangers::cracker::faceBreaker() const
{
    return lawPtr_();
}


Foam::faceBreakerLaw&
Foam::fvMeshTopoChangers::cracker::faceBreaker()
{
    return lawPtr_();
}


void Foam::fvMeshTopoChangers::cracker::topoChange
(
    const polyTopoChangeMap& map
)
{}


void Foam::fvMeshTopoChangers::cracker::mapMesh(const polyMeshMap& map)
{}


void Foam::fvMeshTopoChangers::cracker::distribute
(
    const polyDistributionMap& map
)
{}


// ************************************************************************* //

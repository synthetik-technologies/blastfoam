/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020 Synthetik Applied Technologies: |    Modified original
                            dynamicRefineBalanceBlastFvMesh class
                            to be more appilcable to compressible flows.
                            Improved compatibility with snappyHexMesh.
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

// THIS IS BAD
// The volumes need to be cleared when balancing, but there is no
// public or protected way to do this without using clearOut() which
// will the pointMesh if it exists, so pointFields become invalid.
// This is a work around
#define curTimeIndex_ curTimeIndex_; public:
#include "fvMesh.H"
#undef curTimeIndex_

#include "fvMeshRefiner.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "polyTopoChange.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "pointFields.H"
#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "pointMesh.H"
#include "cellSet.H"
#include "wedgePolyPatch.H"
#include "hexRef3D.H"
#include "RefineBalanceMeshObject.H"
#include "parcelCloud.H"
#include "extrapolatedCalculatedFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshRefiner, 0);
    defineRunTimeSelectionTable(fvMeshRefiner, fvMesh);
    defineRunTimeSelectionTable(fvMeshRefiner, dictionary);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::fvMeshRefiner::count
(
    const PackedBoolList& l,
    const unsigned int val
)
{
    label n = 0;
    forAll(l, i)
    {
        if (l.get(i) == val)
        {
            n++;
        }

        // debug also serves to get-around Clang compiler trying to optimsie
        // out this forAll loop under O3 optimisation
        if (debug)
        {
            Info<< "n=" << n << endl;
        }
    }

    return n;
}


Foam::scalarField
Foam::fvMeshRefiner::maxPointField(const scalarField& pFld) const
{
    scalarField vFld(mesh_.nCells(), -GREAT);

    forAll(mesh_.pointCells(), pointi)
    {
        const labelList& pCells = mesh_.pointCells()[pointi];

        forAll(pCells, i)
        {
            vFld[pCells[i]] = max(vFld[pCells[i]], pFld[pointi]);
        }
    }
    return vFld;
}


Foam::scalarField
Foam::fvMeshRefiner::maxCellField(const scalarField& vFld) const
{
    scalarField pFld(mesh_.nPoints(), -GREAT);

    forAll(mesh_.pointCells(), pointi)
    {
        const labelList& pCells = mesh_.pointCells()[pointi];

        forAll(pCells, i)
        {
            pFld[pointi] = max(pFld[pointi], vFld[pCells[i]]);
        }
    }
    return pFld;
}


Foam::scalarField
Foam::fvMeshRefiner::cellToPoint(const scalarField& vFld) const
{
    scalarField pFld(mesh_.nPoints());

    forAll(mesh_.pointCells(), pointi)
    {
        const labelList& pCells = mesh_.pointCells()[pointi];

        scalar sum = 0.0;
        forAll(pCells, i)
        {
            sum += vFld[pCells[i]];
        }
        pFld[pointi] = sum/pCells.size();
    }
    return pFld;
}


Foam::scalarField Foam::fvMeshRefiner::error
(
    const scalarField& fld,
    const scalar minLevel,
    const scalar maxLevel
) const
{
    scalarField c(fld.size(), -1);

    forAll(fld, i)
    {
        scalar err = min(fld[i]-minLevel, maxLevel-fld[i]);

        if (err >= 0)
        {
            c[i] = err;
        }
    }
    return c;
}


void Foam::fvMeshRefiner::selectRefineCandidates
(
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalarField& vFld,
    PackedBoolList& candidateCell
) const
{
    forAll(vFld, cellI)
    {
        if ((vFld[cellI] >= lowerRefineLevel) && (vFld[cellI] <= upperRefineLevel))
        {
            candidateCell.set(cellI, 1);
        }
    }
}


Foam::labelList Foam::fvMeshRefiner::selectRefineCells
(
    const label maxCells,
    const labelList& maxRefinement,
    const PackedBoolList& candidateCell
) const
{

    // Count current selection
    label nLocalCandidates = count(candidateCell, 1);

    // Collect all cells
    DynamicList<label> candidates(nLocalCandidates);

    forAll(candidateCell, celli)
    {
        if
        (
            cellLevel()[celli] < maxRefinement[celli]
         && candidateCell.get(celli)
        )
        {
            candidates.append(celli);
        }
    }

    candidates.shrink();

    return move(candidates);
}


void Foam::fvMeshRefiner::setMaxCellLevel(labelList& maxCellLevel) const
{
    if (!maxCellLevel.size())
    {
        maxCellLevel.setSize
        (
            mesh_.nCells(),
            dict_.lookup<label>("maxRefinement")
        );
    }

    if (gMin(maxCellLevel) < 0)
    {
        FatalErrorInFunction
            << "Illegal maximum refinement level " << gMin(maxCellLevel) << nl
            << "The maxRefinement should be > 0." << nl
            << exit(FatalError);
    }
    else if (maxCellLevel.size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "Inconsistent number of cells and size of maxCellsLevel "
            << endl
            << abort(FatalError);
    }
}


bool Foam::fvMeshRefiner::preUpdate()
{
    if (canRefine() || canUnrefine())
    {
        HashTable<parcelCloud*> clouds
        (
            mesh_.lookupClass<parcelCloud>()
        );
        forAllIter(HashTable<parcelCloud*>, clouds, iter)
        {
            iter()->storeGlobalPositions();
        }
        return true;
    }
    return false;
}


bool Foam::fvMeshRefiner::canRefine(const bool incr) const
{
    if (!refine_)
    {
        return false;
    }

    const Time& t = mesh_.time();
    if (force_)
    {}
    else if
    (
        t.timeIndex() <= 0
     || t.value() < beginRefine_
     || t.value() > endRefine_
    )
    {
        return false;
    }
    else if ((t.timeIndex() % refineInterval_) > 0)
    {
        return false;
    }

    if (incr)
    {
        nRefinementIterations_++;
    }
    return mesh_.globalData().nTotalCells() < maxCells_;
}


bool Foam::fvMeshRefiner::canUnrefine(const bool incr) const
{
    if (!unrefine_)
    {
        return false;
    }

    const Time& t = mesh_.time();
    if (force_)
    {}
    else if
    (
        t.timeIndex() <= 0
     || t.value() < beginUnrefine_
     || t.value() > endUnrefine_
    )
    {
        return false;
    }
    else if ((t.timeIndex() % unrefineInterval_) > 0)
    {
        return false;
    }

    if (incr)
    {
        nUnrefinementIterations_++;
    }
    return true;
}


bool Foam::fvMeshRefiner::canBalance(const bool incr) const
{
    if (!balancer_.balance())
    {
        return false;
    }

    const Time& t = mesh_.time();

    if (force_)
    {}
    else if
    (
        nRefinementIterations_ <= 0
     || t.value() < beginBalance_
     || t.value() > endBalance_
    )
    {
        return false;
    }
    else if
    (
        (
            max(nRefinementIterations_, nUnrefinementIterations_)
          % balanceInterval_
        ) > 0
    )
    {
        return false;
    }

    // only check if the mesh is unbalanced if everything else is ok
    if (incr)
    {
        nBalanceIterations_++;
    }
    return balancer_.canBalance();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshRefiner::fvMeshRefiner(fvMesh& mesh)
:
    FvMeshRefiner
    (
        mesh,
        IOobject
        (
            typeName,
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    ),

    mesh_(mesh),
    dict_(),

    balancer_(mesh_),

    nRefinementIterations_(0),
    nUnrefinementIterations_(0),
    nBalanceIterations_(0),

    force_(false),
    refine_(true),
    unrefine_(true),

    refineInterval_(1),
    unrefineInterval_(1),
    balanceInterval_(1),

    beginRefine_(0),
    beginUnrefine_(0),
    beginBalance_(0),

    endRefine_(great),
    endUnrefine_(great),
    endBalance_(great),

    maxCells_(labelMax),

    nRefinementBufferLayers_(0),
    nUnrefinementBufferLayers_(0),

    protectedPatches_(),

    dumpLevel_(false),

    isRefining_(false),
    isUnrefining_(false),
    isBalancing_(false),

    V0OldPtr_(nullptr),
    V00OldPtr_(nullptr)
{}


Foam::fvMeshRefiner::fvMeshRefiner
(
    fvMesh& mesh,
    const dictionary& dict,
    const bool force,
    const bool read
)
:
    FvMeshRefiner
    (
        mesh,
        IOobject
        (
            typeName,
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    ),

    mesh_(mesh),
    dict_(dict),

    balancer_
    (
        mesh_,
        dict_.optionalSubDict("loadBalance")
    ),

    nRefinementIterations_(0),
    nUnrefinementIterations_(0),
    nBalanceIterations_(0),

    force_(false),
    refine_(true),
    unrefine_(true),

    refineInterval_(1),
    unrefineInterval_(1),
    balanceInterval_(1),

    beginRefine_(0),
    beginUnrefine_(0),
    beginBalance_(0),

    endRefine_(great),
    endUnrefine_(great),
    endBalance_(great),

    maxCells_(labelMax),

    nRefinementBufferLayers_(0),
    nUnrefinementBufferLayers_(0),

    protectedPatches_(),

    dumpLevel_(false),

    isRefining_(false),
    isUnrefining_(false),
    isBalancing_(false),

    V0OldPtr_(nullptr),
    V00OldPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshRefiner::~fvMeshRefiner()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshRefiner::readDict(const dictionary& dict)
{
    dict_ = dict;

    maxCells_ = dict_.lookupOrDefault("maxCells", labelMax);
    if (maxCells_ <= 0)
    {
        FatalErrorInFunction
            << "Illegal maximum number of cells " << maxCells_ << nl
            << "The maxCells should be > 0." << nl
            << exit(FatalError);
    }

    if (dict.found("nRefinementBufferLayers"))
    {
        nRefinementBufferLayers_ = dict_.lookup<label>("nRefinementBufferLayers");
    }
    else if (dict.found("nBufferLayers"))
    {
        nRefinementBufferLayers_ = dict_.lookup<label>("nBufferLayers");
    }

    if (dict.found("nUnrefinementBufferLayers"))
    {
        nUnrefinementBufferLayers_ =
            dict_.lookup<label>("nUnrefinementBufferLayers");
    }
    else if (dict.found("nBufferLayers"))
    {
        nUnrefinementBufferLayers_ = dict_.lookup<label>("nBufferLayers");
    }

    dumpLevel_ = dict_.lookupOrDefault<bool>("dumpLevel", false);
    protectedPatches_ = dict_.lookupOrDefault("protectedPatches", wordList());

    refine_ = dict_.lookupOrDefault("refine", true);
    unrefine_ = dict_.lookupOrDefault("unrefine", true);

    if (force_)
    {
        refineInterval_ = 1;
        unrefineInterval_ = 1;
        beginUnrefine_ = -great;
        beginBalance_ = -great;
    }
    else
    {
        refineInterval_ = dict_.lookupOrDefault<label>("refineInterval", 1);
        if (refineInterval_ < 0)
        {
            FatalErrorInFunction
                << "Illegal refineInterval " << refineInterval_ << nl
                << "The refineInterval should be >= 1." << nl
                << exit(FatalError);
        }

        unrefineInterval_ =
            dict_.lookupOrDefault<label>("unrefineInterval", refineInterval_);
        if (unrefineInterval_ < 0)
        {
            FatalErrorInFunction
                << "Illegal unrefineInterval " << unrefineInterval_ << nl
                << "The unrefineInterval should be >= 1." << nl
                << exit(FatalError);
        }

        balanceInterval_ =
            dict_.lookupOrDefault<label>("balanceInterval", refineInterval_);
        if (balanceInterval_ < 0)
        {
            FatalErrorInFunction
                << "Illegal balanceInterval " << balanceInterval_ << nl
                << "The balanceInterval should be >= 1." << nl
                << exit(FatalError);
        }

        beginRefine_ = dict_.lookupOrDefault<scalar>("beginRefine", 0.0);
        beginUnrefine_ = dict_.lookupOrDefault<scalar>("beginUnrefine", 0.0);
        beginBalance_ = dict_.lookupOrDefault<scalar>("beginBalance", 0.0);

        endRefine_ = dict_.lookupOrDefault<scalar>("endRefine", great);
        endUnrefine_ = dict_.lookupOrDefault<scalar>("endUnrefine", great);
        endBalance_ = dict_.lookupOrDefault<scalar>("endBalance", great);
    }
}


bool Foam::fvMeshRefiner::balance()
{
    //Part 1 - Call normal update from dynamicRefineBlastFvMesh
    const dictionary& balanceDict(dict_.optionalSubDict("loadBalance"));
    balancer_.read(balanceDict);


    // Part 2 - Load Balancing
    if (canBalance(true))
    {
        isBalancing_ = true;

        //- Save the old volumes so it will be distributed and
        //  resized
        //  We cheat because so we can check which fields
        //  actually need to be mapped
        if (mesh_.V0Ptr_)
        {
            V0OldPtr_ = mesh_.V0Ptr_;
            mesh_.V0Ptr_ = nullptr;
        }
        if (mesh_.V00Ptr_)
        {
            V00OldPtr_ = mesh_.V00Ptr_;
            mesh_.V00Ptr_ = nullptr;
        }

        //- Only clear old volumes if balancing is occurring
        //- Clear V, V0, and V00 since they are not
        //  registered, and therefore are not resized and the
        //  normal mapping does not work.
        //  Instead we save V0/V00 and reset it.

        // The actual fix to this is in progress

        //  THIS IS A PRIVATE FUNCTION OF fvMesh,
        //  but we use a MACRO hack to make it accessible
        mesh_.clearGeom();

        Info<< "Mapping the fields ..." << endl;
        autoPtr<mapDistributePolyMesh> map = balancer_.distribute();

        //- Distribute other data
        distribute(map());

        isBalancing_ = false;

        return true;
    }

    return false;
}


void Foam::fvMeshRefiner::updateMesh(const mapPolyMesh& mpm)
{
    if
    (
        mesh_.foundObject<volScalarField::Internal>("V0_Old")
     || mesh_.foundObject<volScalarField::Internal>("V00_Old")
    )
    {
        //- Only clear old volumes if balancing is occurring
        //- Clear V, V0, and V00 since they are not
        //  registered, and therefore are not resized and the
        //  normal mapping does not work.
        //  Instead we save V0/V00 and reset it.

        // The actual fix to this is in progress

        //  THIS IS A PRIVATE FUNCTION OF fvMesh,
        //  but we use a MACRO hack to make it accessible
        mesh_.clearGeom();
    }
    else
    {
        mesh_.clearGeomNotOldVol();
    }
}


void Foam::fvMeshRefiner::distribute
(
    const mapDistributePolyMesh& map
)
{
    Info<< "Distribute the map ..." << endl;

    //- The volume has been updated, so now we copy back
    //  This also calls V() which will construct the volume
    //  field.
    //  Again, we cheat to access the volume field pointers
    //  This is necessary because the V0 and V00 fields are
    //  not created until the time has advanced and asking for
    //  thermo though V0() or V00() results in a fatal error

    if (V0OldPtr_)
    {
        map.distributeCellData(*V0OldPtr_);
        if (mesh_.V0Ptr_)
        {
            deleteDemandDrivenData(mesh_.V0Ptr_);
        }
        mesh_.V0Ptr_ = V0OldPtr_;
        V0OldPtr_ = nullptr;
    }
    if (V00OldPtr_)
    {
        map.distributeCellData(*V00OldPtr_);
        if (mesh_.V00Ptr_)
        {
            deleteDemandDrivenData(mesh_.V0Ptr_);
        }
        mesh_.V00Ptr_ = V00OldPtr_;
        V00OldPtr_ = nullptr;
    }
}


void Foam::fvMeshRefiner::extendMaxCellLevel
(
    const polyMesh& mesh,
    labelList& cells,
    labelList& maxCellLevel,
    const label level
)
{
    // Mark faces using any marked cell
    boolList markedFace(mesh.nFaces(), false);
    PackedBoolList markedCell(mesh.nCells(), false);

    forAll(cells, i)
    {
        label celli = cells[i];
        markedCell.set(celli, true);
        const cell& cFaces = mesh.cells()[celli];

        forAll(cFaces, i)
        {
            markedFace[cFaces[i]] = true;
        }
    }

    syncTools::syncFaceList(mesh, markedFace, orEqOp<bool>());

    // Update cells using any markedFace
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        if (markedFace[facei])
        {
            markedCell.set(mesh.faceOwner()[facei], 1);
            markedCell.set(mesh.faceNeighbour()[facei], 1);
        }
    }
    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        if (markedFace[facei])
        {
            markedCell.set(mesh.faceOwner()[facei], 1);
        }
    }

    cells.resize(mesh.nCells());
    label i = 0;
    forAll(markedCell, celli)
    {
        if (markedCell.get(celli))
        {
            cells[i++] = celli;
            maxCellLevel[celli] = max(maxCellLevel[celli], level);
        }
    }
    cells.resize(i);
}


void Foam::fvMeshRefiner::extendMarkedCells
(
    PackedBoolList& markedCells,
    const labelList& maxCellLevel,
    const bool isTop,
    const bool force
)
{
    // Mark faces using any marked cell
    boolList markedFace(mesh_.nFaces(), false);

    if (force)
    {
        forAll(markedCells, celli)
        {
            if
            (
                markedCells.get(celli)
             && (maxCellLevel[celli] > cellLevel()[celli] || !isTop)
            )
            {
                const cell& cFaces = mesh_.cells()[celli];

                forAll(cFaces, i)
                {
                    markedFace[cFaces[i]] = true;
                }
            }
        }
    }
    else
    {
        forAll(markedCells, celli)
        {
            if (markedCells.get(celli))
            {
                const cell& cFaces = mesh_.cells()[celli];

                forAll(cFaces, i)
                {
                    markedFace[cFaces[i]] = true;
                }
            }
        }
    }

    syncTools::syncFaceList(mesh_, markedFace, orEqOp<bool>());

    // Update cells using any markedFace
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (markedFace[facei])
        {
            markedCells.set(mesh_.faceOwner()[facei], 1);
            markedCells.set(mesh_.faceNeighbour()[facei], 1);
        }
    }
    for
    (
        label facei = mesh_.nInternalFaces();
        facei < mesh_.nFaces();
        facei++
    )
    {
        if (markedFace[facei])
        {
            markedCells.set(mesh_.faceOwner()[facei], 1);
        }
    }
}


void Foam::fvMeshRefiner::extendMarkedCellsAcrossFaces
(
    PackedBoolList& markedCells
)
{
    // Mark all faces for all marked cells
    const label nFaces = mesh_.nFaces();
    boolList markedFace(nFaces, false);

    // Get mesh cells
    const cellList& meshCells = mesh_.cells();

    // Loop through all cells
    forAll (markedCells, cellI)
    {
        if (markedCells[cellI])
        {
            // This cell is marked, get its faces
            const cell& cFaces = meshCells[cellI];

            forAll (cFaces, i)
            {
                markedFace[cFaces[i]] = true;
            }
        }
    }

    // Snyc the face list across processor boundaries
    syncTools::syncFaceList(mesh_, markedFace, orEqOp<bool>());

    // Get necessary mesh data
    const label nInternalFaces = mesh_.nInternalFaces();
    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    // Internal faces
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        if (markedFace[faceI])
        {
            // Face is marked, mark both owner and neighbour
            const label& own = owner[faceI];
            const label& nei = neighbour[faceI];

            // Mark owner and neighbour cells
            markedCells.set(own, true);
            markedCells.set(nei, true);
        }
    }

    // Boundary faces
    for (label faceI = nInternalFaces; faceI < nFaces; ++faceI)
    {
        if (markedFace[faceI])
        {
            // Face is marked, mark owner
            const label& own = owner[faceI];

            // Mark owner
            markedCells.set(own);
        }
    }
}


void Foam::fvMeshRefiner::extendMarkedCellsAcrossPoints
(
    PackedBoolList& markedCells
)
{
    // Mark all points for all marked cells
    const label nPoints = mesh_.nPoints();
    boolList markedPoint(nPoints, false);

    // Get cell points
    const labelListList& meshCellPoints = mesh_.cellPoints();

    // Loop through all cells
    forAll (markedCells, cellI)
    {
        if (markedCells.get(cellI))
        {
            // This cell is marked, get its points
            const labelList& cPoints = meshCellPoints[cellI];

            forAll (cPoints, i)
            {
                markedPoint[cPoints[i]] = true;
            }
        }
    }

    // Snyc point list across processor boundaries
    syncTools::syncPointList
    (
        mesh_,
        markedPoint,
        orEqOp<bool>(),
        true // Default value
    );

    // Get point cells
    const labelListList& meshPointCells = mesh_.pointCells();

    // Loop through all points
    forAll (markedPoint, pointI)
    {
        if (markedPoint[pointI])
        {
            // This point is marked, mark all of its cells
            const labelList& pCells = meshPointCells[pointI];

            forAll (pCells, i)
            {
                markedCells.set(pCells[i], true);
            }
        }
    }
}


bool Foam::fvMeshRefiner::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    if (dumpLevel_)
    {
        volScalarField scalarCellLevel
        (
            volScalarField::New
            (
                "cellLevel",
                mesh_,
                dimensionedScalar(dimless, 0),
                extrapolatedCalculatedFvPatchField<scalar>::typeName
            )
        );
        forAll(cellLevel(), celli)
        {
            scalarCellLevel[celli] = cellLevel()[celli];
        }
        scalarCellLevel.correctBoundaryConditions();

        pointScalarField scalarPointLevel
        (
            pointScalarField::New
            (
                "pointLevel",
                pointMesh::New(mesh_),
                dimensionedScalar(dimless, 0.0)
            )
        );

        scalarField& sPointLevel = scalarPointLevel.primitiveFieldRef();
        forAll(sPointLevel, pointi)
        {
            sPointLevel[pointi] = pointLevel()[pointi];
        }

        return
            scalarCellLevel.write()
         && scalarPointLevel.write();
    }
    return true;
}

// ************************************************************************* //

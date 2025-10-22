/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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


#include "polyMeshRefiner.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "pointMesh.H"
#include "cellSet.H"
#include "wedgePolyPatch.H"
#include "hexRef3D.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "RefineBalanceMeshObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyMeshRefiner, 0);
    defineRunTimeSelectionTable(polyMeshRefiner, polyMesh);
    defineRunTimeSelectionTable(polyMeshRefiner, dictionary);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::polyMeshRefiner::count
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
Foam::polyMeshRefiner::maxPointField(const scalarField& pFld) const
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
Foam::polyMeshRefiner::maxCellField(const scalarField& vFld) const
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
Foam::polyMeshRefiner::cellToPoint(const scalarField& vFld) const
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


Foam::scalarField Foam::polyMeshRefiner::error
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


void Foam::polyMeshRefiner::selectRefineCandidates
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


Foam::labelList Foam::polyMeshRefiner::selectRefineCells
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

    return candidates;
}


void Foam::polyMeshRefiner::setMaxCellLevel(labelList& maxCellLevel) const
{
    if (!returnReduce(maxCellLevel.size(), sumOp<label>()))
    {
        if (!dict_.found("maxRefinement"))
        {
            FatalIOErrorInFunction(dict_)
                << "maxRefinement was not specified" << endl
                << abort(FatalError);
        }
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


bool Foam::polyMeshRefiner::preUpdate()
{
    if (canRefine() || canUnrefine())
    {
        return true;
    }
    return false;
}


bool Foam::polyMeshRefiner::canRefine(const bool incr) const
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


bool Foam::polyMeshRefiner::canUnrefine(const bool incr) const
{
    // Make sure location mapper is cleared since removed points do not
    // need to be updated
    locationMapper::NewRef(mesh_).clearOut();

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyMeshRefiner::polyMeshRefiner(polyMesh& mesh)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            mesh.facesInstance(),
            mesh
        )
    ),

    mesh_(mesh),
    dict_(),

    nRefinementIterations_(0),
    nUnrefinementIterations_(0),

    force_(false),
    refine_(true),
    unrefine_(true),

    refineInterval_(1),
    unrefineInterval_(1),

    beginRefine_(0),
    beginUnrefine_(0),

    endRefine_(great),
    endUnrefine_(great),

    maxCells_(labelMax),

    nRefinementBufferLayers_(0),
    nUnrefinementBufferLayers_(0),

    protectedPatches_(),

    isRefining_(false),
    isUnrefining_(false)
{
    locationMapper::New(mesh_);
}


Foam::polyMeshRefiner::polyMeshRefiner
(
    polyMesh& mesh,
    const dictionary& dict,
    const bool read
)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),
    dict_(dict),

    nRefinementIterations_(0),
    nUnrefinementIterations_(0),

    force_(dict.lookupOrDefault("force", false)),
    refine_(true),
    unrefine_(true),

    refineInterval_(1),
    unrefineInterval_(1),

    beginRefine_(0),
    beginUnrefine_(0),

    endRefine_(great),
    endUnrefine_(great),

    maxCells_(labelMax),

    nRefinementBufferLayers_(0),
    nUnrefinementBufferLayers_(0),

    protectedPatches_(),

    isRefining_(false),
    isUnrefining_(false)
{
    locationMapper::New(mesh_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyMeshRefiner::~polyMeshRefiner()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyMeshRefiner::readDict(const dictionary& dict)
{
    if (&dict_ != &dict)
    {
        dict_ = dict;
    }

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

    protectedPatches_ = dict_.lookupOrDefault("protectedPatches", wordList());

    refine_ = dict_.lookupOrDefault("refine", true);
    unrefine_ = dict_.lookupOrDefault("unrefine", true);

    if (force_)
    {
        refineInterval_ = 1;
        unrefineInterval_ = 1;
        beginUnrefine_ = -great;
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

        beginRefine_ = dict_.lookupOrDefault<scalar>("beginRefine", 0.0);
        beginUnrefine_ = dict_.lookupOrDefault<scalar>("beginUnrefine", 0.0);

        endRefine_ = dict_.lookupOrDefault<scalar>("endRefine", great);
        endUnrefine_ = dict_.lookupOrDefault<scalar>("endUnrefine", great);
    }
}


void Foam::polyMeshRefiner::topoChange(const polyTopoChangeMap& map)
{
    const locationMapper& locMapper = locationMapper::New(mesh_);
    const wordHashSet& interpolatedFields =
        locMapper.interpolatedFields();

    const labelList& pointMap = map.pointMap();

    forAllConstIter(wordHashSet, interpolatedFields, iter)
    {
        const word& fieldName = iter.key();
        if (mesh_.foundObject<pointIOField>(fieldName))
        {
            pointIOField& points =
                mesh_.lookupObjectRef<pointIOField>(fieldName);

            pointField newPoints(pointMap.size(), Zero);
            forAll(newPoints, pointi)
            {
                label oldPointi = pointMap[pointi];

                if (oldPointi >= 0)
                {
                    newPoints[pointi] = points[oldPointi];
                }
            }
            locMapper.interpolateMidPoints(newPoints);
            points.transfer(newPoints);
        }
        else
        {
            WarningInFunction
                << fieldName << " is not a registered pointField. Not mapping"
                << endl;
        }
    }
}


void Foam::polyMeshRefiner::mapMesh(const polyMeshMap&)
{
    NotImplemented;
}


void Foam::polyMeshRefiner::distribute
(
    const polyDistributionMap& map
)
{
    const locationMapper& locMapper = locationMapper::New(mesh_);
    const wordHashSet& interpolatedFields =
        locMapper.interpolatedFields();

    forAllConstIter(wordHashSet, interpolatedFields, iter)
    {
        const word& fieldName = iter.key();
        if (mesh_.foundObject<pointIOField>(fieldName))
        {
             map.distributePointData
            (
                mesh_.lookupObjectRef<pointIOField>(fieldName)
            );
        }
        else
        {
            WarningInFunction
                << fieldName << " is not a registered pointField. Not mapping"
                << endl;
        }
    }
}


void Foam::polyMeshRefiner::extendMaxCellLevel
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


void Foam::polyMeshRefiner::extendMarkedCells
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
            markedCells.set(mesh_.faceOwner()[facei]);
            markedCells.set(mesh_.faceNeighbour()[facei]);
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
            markedCells.set(mesh_.faceOwner()[facei]);
        }
    }
}


void Foam::polyMeshRefiner::extendMarkedCellsAcrossFaces
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
            // Mark owner and neighbour cells
            markedCells.set(owner[faceI]);
            markedCells.set(neighbour[faceI]);
        }
    }

    // Boundary faces
    for (label faceI = nInternalFaces; faceI < nFaces; ++faceI)
    {
        if (markedFace[faceI])
        {
            // Mark owner
            markedCells.set(owner[faceI]);
        }
    }
}


void Foam::polyMeshRefiner::extendMarkedCellsAcrossPoints
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
                markedCells.set(pCells[i]);
            }
        }
    }
}


bool Foam::polyMeshRefiner::write(const bool write) const
{
    return true;
}


bool Foam::polyMeshRefiner::writeVolFields(const fvMesh& mesh) const
{
    bool writeOk = true;
    {
        tmp<volScalarField> tscalarCellLevel
        (
            volScalarField::New
            (
                "cellLevel",
                mesh,
                dimensionedScalar(dimless, 0),
                extrapolatedCalculatedFvPatchField<scalar>::typeName
            )
        );
        tscalarCellLevel.ref().primitiveFieldRef() =
            scalarList(this->cellLevel());
        tscalarCellLevel.ref().correctBoundaryConditions();
        writeOk = tscalarCellLevel().write() && writeOk;
    }

    {
        tmp<pointScalarField> tscalarPointLevel
        (
            pointScalarField::New
            (
                "pointLevel",
                pointMesh::New(mesh),
                dimensionedScalar(dimless, 0.0)
            )
        );
        tscalarPointLevel.ref().primitiveFieldRef() =
            scalarList(this->pointLevel());
        writeOk = tscalarPointLevel().write() && writeOk;
    }
    return writeOk;
}


// ************************************************************************* //

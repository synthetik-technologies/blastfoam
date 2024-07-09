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

bool Foam::fvMeshRefiner::canBalance(const bool incr) const
{
    if (!balancer_.balance())
    {
        return false;
    }

    const Time& t = mesh_.time();

    if (refiner_->force_)
    {}
    else if
    (
        refiner_->nRefinementIterations_ <= 0
     || t.value() < beginBalance_
     || t.value() > endBalance_
    )
    {
        return false;
    }
    else if
    (
        (
            max(refiner_->nRefinementIterations_, refiner_->nUnrefinementIterations_)
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
    return returnReduce(balancer_.canBalance(), orOp<bool>());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshRefiner::fvMeshRefiner
(
    const word& refinerType,
    fvMesh& mesh
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

    refiner_(polyMeshRefiner::New(refinerType, mesh)),
    balancer_(mesh_),

    nBalanceIterations_(0),
    balanceInterval_(1),
    beginBalance_(0),
    endBalance_(great),

    dumpLevel_(false),

    V0OldPtr_(nullptr),
    V00OldPtr_(nullptr)
{}


Foam::fvMeshRefiner::fvMeshRefiner
(
    const word& refinerType,
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
    refiner_(polyMeshRefiner::New(refinerType, mesh, dict, force, read)),

    balancer_
    (
        mesh_,
        refiner_->dict_.optionalSubDict("loadBalance")
    ),

    nBalanceIterations_(0),
    balanceInterval_(1),
    beginBalance_(0),
    endBalance_(great),

    dumpLevel_(false),

    V0OldPtr_(nullptr),
    V00OldPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshRefiner::~fvMeshRefiner()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshRefiner::readDict(const dictionary& dict)
{
    refiner_->readDict(dict);

    dumpLevel_ = refiner_->dict_.lookupOrDefault<bool>("dumpLevel", false);

    if (refiner_->force_)
    {
        beginBalance_ = -great;
    }
    else
    {
        balanceInterval_ =
            refiner_->dict_.lookupOrDefault<label>("balanceInterval", 1);
        if (balanceInterval_ < 0)
        {
            FatalErrorInFunction
                << "Illegal balanceInterval " << balanceInterval_ << nl
                << "The balanceInterval should be >= 1." << nl
                << exit(FatalError);
        }

        beginBalance_ = refiner_->dict_.lookupOrDefault<scalar>("beginBalance", 0.0);
        endBalance_ = refiner_->dict_.lookupOrDefault<scalar>("endBalance", great);
    }
}


bool Foam::fvMeshRefiner::refine
(
    const scalarField& error,
    const labelList& maxCellLevel,
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalar unrefineLevel
)
{
    bool hasChanged =
        refiner_->refine
        (
            error,
            maxCellLevel,
            lowerRefineLevel,
            upperRefineLevel,
            unrefineLevel
        );

    if (balance())
    {
        hasChanged = true;

        mesh_.topoChanging(true);

        // Reset moving flag (if any). If not using inflation we'll not
        // move, if are using inflation any follow on movePoints will set
        // it.
        mesh_.moving(false);

        // Make sure all processors have the correct instance
        mesh_.setInstance(mesh_.time().timeName());
        mesh_.polyMesh::instance() = mesh_.time().timeName();
    }

    return hasChanged;
}


bool Foam::fvMeshRefiner::balance()
{
    //Part 1 - Reread the balance dictionary
    const dictionary& balanceDict(refiner_->dict_.optionalSubDict("loadBalance"));
    balancer_.read(balanceDict);

    // Part 2 - Load Balancing
    if (canBalance(true))
    {
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
        balancer_.distribute();

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
    // else
    // {
    //     mesh_.clearGeomNotOldVol();
    // }

    if (refiner_->isBalancing())
    {
        return;
    }
    const locationMapper& locMapper = locationMapper::New(mesh_);
    const wordHashSet& interpolatedPointFields =
        locMapper.interpolatedPointFields();
    forAllConstIter(wordHashSet, interpolatedPointFields, iter)
    {
        const word& fieldName = iter.key();
        if (mesh_.foundObject<pointVectorField>(fieldName))
        {
            pointVectorField& points =
                mesh_.lookupObjectRef<pointVectorField>(fieldName);
            points.correctBoundaryConditions();
            locMapper.interpolateMidPoints(points.primitiveFieldRef());
            pointConstraints::New(points.mesh()).setPatchFields(points);
        }
        else
        {
            WarningInFunction
                << fieldName << " is not a registered pointVectorField. "
                << "Not mapping" << endl;
        }
    }
}


void Foam::fvMeshRefiner::distribute
(
    const mapDistributePolyMesh& map
)
{
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


bool Foam::fvMeshRefiner::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    bool writeOK = balancer_.write(write) && refiner_->write(write);

    if (dumpLevel_ && write)
    {
        volScalarField scalarCellLevel
        (
            IOobject
            (
                "cellLevel",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless, 0),
            extrapolatedCalculatedFvPatchField<scalar>::typeName
        );
        scalarCellLevel.primitiveFieldRef() = scalarList(refiner_->cellLevel());
        scalarCellLevel.correctBoundaryConditions();

        pointScalarField scalarPointLevel
        (
            IOobject
            (
                "pointLevel",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointMesh::New(mesh_),
            dimensionedScalar(dimless, 0.0)
        );
        scalarPointLevel.primitiveFieldRef() = scalarList(refiner_->pointLevel());

        return
            writeOK
         && scalarCellLevel.write()
         && scalarPointLevel.write();
    }
    return writeOK;
}

// ************************************************************************* //

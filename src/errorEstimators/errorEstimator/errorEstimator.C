/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "errorEstimator.H"
#include "coupledMaxErrorFvPatchScalarField.H"
#include "mappedPatchBase.H"
#include "timeControlFunctionObject.H"
#include "probes.H"
#include "blastProbes.H"
#include "cellSet.H"
#include "meshSizeObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(errorEstimator, 0);
    defineRunTimeSelectionTable(errorEstimator, dictionary);
}


// * * * * * * * * * * * * * Protected Member Function * * * * * * * * * * * //

Foam::volScalarField& Foam::errorEstimator::lookupOrConstructError
(
    const fvMesh& mesh
) const
{
    word errorName(IOobject::groupName("error", name_));

    if (!mesh.foundObject<volScalarField>(errorName))
    {
        wordList boundaryTypes(mesh.boundaryMesh().size(), "zeroGradient");
        forAll(boundaryTypes, patchi)
        {
            if
            (
                isA<mappedPatchBase>(mesh.boundary()[patchi])
            )
            {
                boundaryTypes[patchi] =
                    coupledMaxErrorFvPatchScalarField::typeName;
            }
            if (debug)
            {
                Pout<< "Patch:" << mesh.boundary()[patchi].patch().name() <<nl
                    << " cellType:" << boundaryTypes[patchi] << endl;
            }
        }

        volScalarField* fPtr =
            new volScalarField
            (
                IOobject
                (
                    errorName,
                    mesh.time().name(),
                    mesh
                ),
                mesh,
                0.0,
                boundaryTypes
            );
        fPtr->store(fPtr);
    }
    return mesh.lookupObjectRef<volScalarField>(errorName);
}

bool Foam::errorEstimator::updateCurTimeIndex(const bool unset) const
{
    if (force_)
    {
        return false;
    }
    if (unset)
    {
        curTimeIndex_--;
        return false;
    }
    if (curTimeIndex_ != mesh_.time().timeIndex())
    {
        curTimeIndex_ = mesh_.time().timeIndex();
        return false;
    }
    return true;
}


Foam::labelHashSet Foam::errorEstimator::errorCells() const
{
    if (cZones_.size())
    {
        const cellZoneList& zones = mesh_.cellZones();
        labelHashSet cells;
        if (hasDefault_)
        {
            cellSet allCells
            (
                mesh_,
                "allCells",
                IOobject::NO_READ,
                IOobject::NO_WRITE
            );
            forAll(zones, i)
            {
                allCells.insert(zones[i]);
            }
            allCells.invert(mesh_.nCells());
            cells = allCells;
        }

        forAll(cZones_, i)
        {
            cells.insert(mesh_.cellZones()[cZones_[i]]);
        }
        return cells;
    }
    return identityMap(mesh_.nCells());
}


template<>
bool Foam::errorEstimator::getFieldValueType<Foam::scalar>
(
    const word& name,
    volScalarField& f,
    const labelHashSet& eCells
) const
{
    typedef GeometricField<scalar, fvPatchField, volMesh> thisType;

    if (mesh_.foundObject<thisType>(name))
    {
        const thisType& x = mesh_.lookupObject<thisType>(name);
        forAllConstIter(labelHashSet, eCells, iter)
        {
            const label celli = iter.key();
            f[celli] = x[celli];
        }
        return true;
    }
    return false;
}


void Foam::errorEstimator::readCellZones(const dictionary& dict)
{
    wordList zones;
    if (dict.found("zones"))
    {
        dict.readIfPresent("zones", zones);
    }
    else if (dict.found("zone"))
    {
        zones.setSize(1);
        dict.readIfPresent("zone", zones[0]);
    }
    hasDefault_ = false;
    forAll(zones, i)
    {
        if (zones[i] == "default")
        {
            hasDefault_ = true;
        }
        else if (!mesh_.cellZones().found(zones[i]))
        {
            FatalIOErrorInFunction(dict)
                << zones[i] << " is not a valid cell zone. Use \"default\" "
                << "to use all cell not belonging to a cell zone or " << nl
                << mesh_.cellZones().toc() << endl
                << abort(FatalIOError);
        }
        else
        {
            cZones_.append(zones[i]);
        }
    }
}

void Foam::errorEstimator::readMaxRefinement(const dictionary& dict)
{
    if (dict.found("maxZoneRefinement"))
    {
        dict.lookup("maxZoneRefinement") >> maxLevel_;
        forAllConstIter(HashTable<label>, maxLevel_, iter)
        {
            if (!mesh_.cellZones().found(iter.key()))
            {
                FatalIOErrorInFunction(dict)
                    << iter.key() << " is not a cell zone, valid option are"
                    << nl
                    << mesh_.cellZones().toc() << endl
                    << abort(FatalIOError);
            }
        }
    }
    if (dict.found("minZoneDx"))
    {
        dict.lookup("minZoneDx") >> minDx_;
        forAllConstIter(HashTable<scalar>, minDx_, iter)
        {
            if (!mesh_.cellZones().found(iter.key()))
            {
                FatalIOErrorInFunction(dict)
                    << iter.key() << " is not a cell zone, valid option are"
                    << nl
                    << mesh_.cellZones().toc() << endl
                    << abort(FatalIOError);
            }
        }
    }
    forAllConstIter(HashTable<label>, maxLevel_, iter)
    {
        if (minDx_.found(iter.key()))
        {
            FatalIOErrorInFunction(dict)
                << "Both maxRefinement and minDx were specified for "
                << iter.key() << ", only one should be specified" << endl
                << abort(FatalIOError);
        }
    }

    if (dict.found("maxRefinement"))
    {
        defaultMaxLevel_ = dict.lookup<label>("maxRefinement");
        defaultMinDx_ = -1;
    }
    else if (dict.found("minDx"))
    {
        defaultMaxLevel_ = dict.lookup<scalar>("minDx");
        defaultMinDx_ = -1;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Either maxRefinement or minDx must be specified" << endl
            << abort(FatalIOError);
    }

    override_ = dict.lookupOrDefault("override", false);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorEstimator::errorEstimator
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, name),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            debug ? IOobject::AUTO_WRITE : IOobject::NO_WRITE,
            name == word::null ? true : false
        )
    ),
    mesh_(mesh),
    name_(name),
    error_(lookupOrConstructError(mesh)),
    lowerRefine_(0.0),
    lowerUnrefine_(0.0),
    upperRefine_(0.0),
    upperUnrefine_(0.0),
    maxLevel_(-1),
    minDx_(-1),
    override_(false),
    cZones_(),
    hasDefault_(false),
    refineProbes_(dict.lookupOrDefault("refineProbes", true)),
    force_(false),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimator::~errorEstimator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::errorEstimator::read(const dictionary& dict)
{
    lowerRefine_ = dict.lookup<scalar>("lowerRefineLevel");
    lowerUnrefine_ = dict.lookup<scalar>("unrefineLevel");
    upperRefine_ = dict.lookupOrDefault("upperRefineLevel", great);
    upperUnrefine_ = dict.lookupOrDefault("upperUnrefineLevel", great);

    readCellZones(dict);
    readMaxRefinement(dict);


}


void Foam::errorEstimator::getFieldValue
(
    const word& name,
    volScalarField& f,
    const labelHashSet& eCells
) const
{
    bool found = false;
    found = found || this->getFieldValueType<scalar>(name, f, eCells);
    found = found || this->getFieldValueType<vector>(name, f, eCells);
    found = found || this->getFieldValueType<symmTensor>(name, f, eCells);
    found = found || this->getFieldValueType<sphericalTensor>(name, f, eCells);
    found = found || this->getFieldValueType<tensor>(name, f, eCells);

    if (!found && f.time().timeIndex() > 0)
    {
        FatalErrorInFunction
            << name << " is not a registered field" << endl
            << abort(FatalError);
    }
}


void Foam::errorEstimator::normalize
(
    volScalarField& error,
    const labelHashSet& eCells
)
{
    error_.correctBoundaryConditions();
    forAllConstIter(labelHashSet, eCells, iter)
    {
        const label celli = iter.key();
        if
        (
            error[celli] < lowerUnrefine_
         || error[celli] > upperUnrefine_
        )
        {
            error[celli] = -1.0;
        }
        else if
        (
            error[celli] > lowerRefine_
         && error[celli] < upperRefine_
        )
        {
            error[celli] = 1.0;
        }
        else
        {
            error[celli] = 0.0;
        }
    }

    volScalarField::Boundary& berror = error.boundaryFieldRef();
    forAll(berror, patchi)
    {
        fvPatchScalarField& perror = berror[patchi];
        const labelList& faceCells = perror.patch().faceCells();

        forAll(perror, facei)
        {
            const label celli = faceCells[facei];
            if (eCells.found(celli))
            {
                if
                (
                    perror[facei] < lowerUnrefine_
                 || perror[facei] > upperUnrefine_
                )
                {
                    error[celli] = max(error[celli], -1.0);
                }
                else if
                (
                    perror[facei] > lowerRefine_
                 && perror[facei] < upperRefine_
                )
                {
                    error[celli] = max(error[celli], 1.0);
                }
                else
                {
                    error[celli] = max(error[celli], 0.0);
                }
            }
        }
    }

    if (!refineProbes_)
    {
        return;
    }

    const functionObjectList& funcs(mesh_.time().functionObjects());
    labelList map;
    forAll(funcs, i)
    {
        vectorField pts;
        if (isA<probes>(funcs[i]))
        {
            const probes& p(refCast<const probes>(funcs[i]));
            pts = p;

        }
        if (isA<blastProbes>(funcs[i]))
        {
            const blastProbes& p(refCast<const blastProbes>(funcs[i]));
            pts = p;
        }
        forAll(pts, j)
        {
            const label celli =
                mesh_.findCell(pts[j], polyMesh::FACE_PLANES);
            if (celli >= 0 && eCells.found(celli))
            {
                error[celli] = 1.0;
            }
        }
    }
}

Foam::labelList Foam::errorEstimator::maxRefinement() const
{
    const labelHashSet& eCells = errorCells();
    labelList maxLevel(mesh_.nCells(), 0);
    if (defaultMaxLevel_ >= 0)
    {
        maxLevel = defaultMaxLevel_;
        if (cZones_.size())
        {
            const labelHashSet& eCells = errorCells();
            maxLevel = 0;
            forAllConstIter(labelHashSet, eCells, iter)
            {
                maxLevel[iter.key()] = defaultMaxLevel_;
            }
        }
    }

    forAllConstIter(HashTable<label>, maxLevel_, iter)
    {
        const cellZone& zone = mesh_.cellZones()[iter.key()];
        forAll(zone, ci)
        {
            const label celli = zone[ci];
            maxLevel[celli] = max(maxLevel[celli], iter());
        }
    }

    if (defaultMaxLevel_ < 0 || minDx_.size())
    {
        const labelIOList& cellLevel
        (
            mesh_.lookupObject<labelIOList>("cellLevel")
        );
        const scalarField& dx(meshSizeObject::New(mesh_).dx());
        if (defaultMaxLevel_ < 0)
        {
            forAllConstIter(labelHashSet, eCells, iter)
            {
                const label celli = iter.key();
                label level = cellLevel[celli];
                if (dx[celli] > defaultMinDx_ && error_[celli] > 0)
                {
                    level++;
                }
                maxLevel[celli] = (maxLevel[celli], level);
            }
        }

        forAllConstIter(HashTable<scalar>, minDx_, iter)
        {
            const cellZone& zone = mesh_.cellZones()[iter.key()];
            forAll(zone, ci)
            {
                const label celli = zone[ci];
                label level = cellLevel[celli];
                if (dx[celli] > iter() && error_[celli] > 0)
                {
                    level++;
                }
                maxLevel[celli] = max(maxLevel[celli], level);
            }
        }
    }

    return maxLevel;
}


bool Foam::errorEstimator::writeData(Ostream&) const
{
    if (debug)
    {
        const_cast<errorEstimator&>(*this).update();
        volScalarField maxLevel
        (
            volScalarField::New
            (
                "maxLevel",
                mesh_,
                0.0
            )
        );
        labelList mr(maxRefinement());
        forAll(mr, celli)
        {
            maxLevel[celli] = mr[celli];
        }
        maxLevel.write();

        if (defaultMinDx_ > 0 || minDx_.size())
        {
            const meshSizeObject& mso = meshSizeObject::New(mesh_);
            const_cast<meshSizeObject&>(mso).movePoints();
            mso.dx(mesh_)().write();
            mso.dX(mesh_)().write();
        }

        return error_.write();
    }
    return true;
}

// ************************************************************************* //

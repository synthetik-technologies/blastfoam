/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2022
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
                boundaryTypes[patchi] = coupledMaxErrorFvPatchScalarField::typeName;
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
                    mesh.time().timeName(),
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
            mesh.time().timeName(),
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
    refineProbes_(dict.lookupOrDefault("refineProbes", true)),
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

    if (dict.found("maxRefinement"))
    {
        maxLevel_ = dict.lookup<label>("maxRefinement");
        minDx_ = -1;
    }
    else if (dict.found("minDx"))
    {
        minDx_ = dict.lookup<scalar>("minDx");
        maxLevel_ = -1;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Either maxRefinement or minDx must be specified" << endl
            << abort(FatalIOError);
    }
}


void Foam::errorEstimator::getFieldValue(const word& name, volScalarField& f) const
{
    this->getFieldValueType<scalar>(name, f);
    this->getFieldValueType<vector>(name, f);
    this->getFieldValueType<symmTensor>(name, f);
    this->getFieldValueType<sphericalTensor>(name, f);
    this->getFieldValueType<tensor>(name, f);
}


void Foam::errorEstimator::normalize(volScalarField& error)
{
    forAll(error, celli)
    {
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
            label celli = mesh_.findCell(pts[j], polyMesh::FACE_PLANES);
            if (celli >= 0)
            {
                error[celli] = 1.0;
            }
        }
    }
}

Foam::labelList Foam::errorEstimator::maxRefinement() const
{
    if (maxLevel_ >= 0 || !mesh_.foundObject<labelIOList>("cellLevel"))
    {
        return labelList(mesh_.nCells(), maxLevel_);
    }

    const labelIOList& cellLevel
    (
        mesh_.lookupObject<labelIOList>("cellLevel")
    );
    labelList maxLevel(cellLevel);
    vector validD(mesh_.geometricD());
    for (label cmpti = 0; cmpti < 3; cmpti++)
    {
        if (validD[cmpti] < 0)
        {
            validD[cmpti] = great;
        }
    }
    const scalarField& dx(meshSizeObject::New(mesh_).dx());

    forAll(dx, celli)
    {
        if (dx[celli] > minDx_ && error_[celli] > 0)
        {
            maxLevel[celli]++;
        }
    };
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

        if (minDx_ > 0)
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

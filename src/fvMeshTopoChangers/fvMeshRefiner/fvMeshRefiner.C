/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022-2025
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

#define list fvMeshRefiner
#include "fvMesh.H"

#include "fvMeshRefiner.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "zeroGradientFvPatchField.H"
#include "syncTools.H"
#include "pointFields.H"
#include "pointMesh.H"
#include "fvcPointInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(fvMeshRefiner, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, fvMeshRefiner, fvMesh);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::fvMeshRefiner::fvMeshRefiner
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    fvMeshRefiner(mesh, dict, dict.lookup<word>("refiner"))
{}


Foam::fvMeshTopoChangers::fvMeshRefiner::fvMeshRefiner
(
    fvMesh& mesh,
    const dictionary& dict,
    const word& refinerType
)
:
    fvMeshTopoChanger(mesh),

    error_(errorEstimator::New(mesh, dict)),

    refiner_(polyMeshRefiner::New(refinerType, mesh, dict, true)),

    dumpLevel_(dict.lookupOrDefault("dumpLevel", true))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::fvMeshRefiner::~fvMeshRefiner()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::fvMeshRefiner::update()
{
    error_->update();

    bool updated = refine
    (
        error_->error(),
        error_->maxRefinement(),
        sqrt(small),
        great,
        -sqrt(small)
    );
    if (updated)
    {
        mesh().moving_ = false;
    }
    return updated;
}


bool Foam::fvMeshTopoChangers::fvMeshRefiner::refine
(
    const scalarField& error,
    const labelList& maxCellLevel,
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalar unrefineLevel
)
{
    bool updated = refiner_->refine
    (
        error,
        maxCellLevel,
        lowerRefineLevel,
        upperRefineLevel,
        unrefineLevel
    );
    if (updated)
    {
        mesh().moving_ = false;
    }
    return updated;
}


void Foam::fvMeshTopoChangers::fvMeshRefiner::topoChange(const polyTopoChangeMap& map)
{
    refiner_->topoChange(map);
}


void Foam::fvMeshTopoChangers::fvMeshRefiner::mapMesh(const polyMeshMap& map)
{
    refiner_->mapMesh(map);
}


void Foam::fvMeshTopoChangers::fvMeshRefiner::distribute(const polyDistributionMap& map)
{
    refiner_->distribute(map);
}


bool Foam::fvMeshTopoChangers::fvMeshRefiner::write(const bool write) const
{
    if (dumpLevel_ && write)
    {
        volScalarField scalarCellLevel
        (
            volScalarField::New
            (
                "cellLevel",
                mesh(),
                dimensionedScalar(dimless, 0),
                zeroGradientFvPatchField<scalar>::typeName
            )
        );
        scalarCellLevel.primitiveFieldRef() = scalarList(refiner_->cellLevel());
        scalarCellLevel.correctBoundaryConditions();

        pointScalarField scalarPointLevel
        (
            pointScalarField::New
            (
                "pointLevel",
                pointMesh::New(mesh()),
                dimensionedScalar(dimless, 0.0)
            )
        );
        scalarPointLevel.primitiveFieldRef() = scalarList(refiner_->pointLevel());
        return scalarCellLevel.write() && scalarPointLevel.write();
    }
    return true;
}

// ************************************************************************* //

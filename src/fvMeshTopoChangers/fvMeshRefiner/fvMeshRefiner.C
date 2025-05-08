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
    defineTypeNameAndDebug(fvMeshRefiner, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, fvMeshRefiner, fvMesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshRefiner::fvMeshRefiner
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    fvMeshRefiner(mesh, dict, dict.lookup<word>("refiner"))
{}


Foam::fvMeshRefiner::fvMeshRefiner
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

Foam::fvMeshRefiner::~fvMeshRefiner()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshRefiner::update()
{
    error_->update();
    error_->error().correctBoundaryConditions();
    return refine
    (
        error_->error(),
        error_->maxRefinement(),
        sqrt(small),
        great,
        -sqrt(small)
    );
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
    return refiner_->refine
    (
        error,
        maxCellLevel,
        lowerRefineLevel,
        upperRefineLevel,
        unrefineLevel
    );
}


void Foam::fvMeshRefiner::topoChange(const polyTopoChangeMap& map)
{
    refiner_->topoChange(map);
}


void Foam::fvMeshRefiner::mapMesh(const polyMeshMap& map)
{
    refiner_->mapMesh(map);
}


void Foam::fvMeshRefiner::distribute(const polyDistributionMap& map)
{
    refiner_->distribute(map);
}


bool Foam::fvMeshRefiner::write(const bool write) const
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

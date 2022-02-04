/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "burstFvPatchBase.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::burstFvPatchBase::updateDeltas()
{
    const fvMesh& mesh(this->patch_.boundaryMesh().mesh());
    surfaceScalarField& weights =
        const_cast<surfaceScalarField&>(mesh.weights());
    surfaceScalarField& deltaCoeffs =
        const_cast<surfaceScalarField&>(mesh.deltaCoeffs());
    surfaceScalarField& nonOrthDeltaCoeffs =
        const_cast<surfaceScalarField&>(mesh.nonOrthDeltaCoeffs());
    surfaceVectorField& nonOrthCorrectionVectors =
        const_cast<surfaceVectorField&>(mesh.nonOrthCorrectionVectors());

    const label patchI = patch_.index();
    const burstFvPatchBase& burstPatch
    (
        dynamicCast<const burstFvPatchBase>(patch_)
    );
    const vectorField pdelta(patch_.fvPatch::delta());
    const scalarField pmagDelta(mag(pdelta));
    const vectorField n
    (
        mesh.Sf().boundaryField()[patchI]
        /mesh.magSf().boundaryField()[patchI]
    );
    scalarField& pweights = weights.boundaryFieldRef()[patchI];
    scalarField& pdeltaCoeffs = deltaCoeffs.boundaryFieldRef()[patchI];
    scalarField& pnonOrthDeltaCoeffs =
        nonOrthDeltaCoeffs.boundaryFieldRef()[patchI];
    vectorField& pnonOrthCorrectionVectors =
        nonOrthCorrectionVectors.boundaryFieldRef()[patchI];

    burstPatch.makeWeights(pweights);
    pweights = intact_ + (1.0 - intact_)*pweights;
    pdeltaCoeffs = 1.0/pmagDelta;
    pnonOrthDeltaCoeffs = 1.0/max(patch_.nf() & pdelta, 0.05*pmagDelta);
    pnonOrthCorrectionVectors = n - pdelta*pnonOrthDeltaCoeffs;
}


void Foam::burstFvPatchBase::setIntact(const scalarField& intact)
{
    burstPolyPatch_.setIntact(intact);
    // Do not update deltas since the mesh may not be constructed yet
}


void Foam::burstFvPatchBase::update
(
    const scalarField& p,
    const scalarField& impulse
)
{
    if (burstPolyPatch_.update(p, impulse))
    {
        updateDeltas();
    }
}
// ************************************************************************* //

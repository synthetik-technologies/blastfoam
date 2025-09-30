/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2020
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

#include "densityGradientErrorEstimator.H"
#include "fvc.H"
#include "meshSizeObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace errorEstimators
{
    defineTypeNameAndDebug(gradient, 0);
    addToRunTimeSelectionTable(errorEstimator, gradient, dictionary);

    defineTypeNameAndDebug(densityGradient, 0);
    addToRunTimeSelectionTable(errorEstimator, densityGradient, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorEstimators::gradient::gradient
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
:
    errorEstimator(mesh, dict, name),
    fieldName_(this->lookupFieldName(dict, typeName))
{
    this->read(dict);
}


Foam::errorEstimators::gradient::gradient
(
    const word& fieldName,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
:
    errorEstimator(mesh, dict, name),
    fieldName_(fieldName)
{
    this->read(dict);
}



Foam::errorEstimators::densityGradient::densityGradient
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
:
    gradient("rho", mesh, dict, name)
{}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimators::gradient::~gradient()
{}


Foam::errorEstimators::densityGradient::~densityGradient()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::errorEstimators::gradient::update(const bool scale)
{
    if (updateCurTimeIndex(!scale))
    {
        return;
    }

    const labelHashSet& eCells = this->errorCells();

    tmp<volScalarField> tx(this->getFieldValue(fieldName_, eCells));
    const volScalarField& x = tx();

    tmp<volVectorField> tgradX(fvc::grad(x));
    const volVectorField& gradX = tgradX();

    const scalarField& dL(meshSizeObject::New(mesh_).dx());

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const label nInternalFaces = mesh_.nInternalFaces();
    error_ = 0.0;

    const vector solutionD((vector(mesh_.geometricD()) + vector::one)/2.0);

    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        const vector dr = mesh_.C()[nei] - mesh_.C()[own];
        const scalar magdr = mag(dr);

        const bool foundOwn = eCells.found(own);
        const bool foundNei= eCells.found(nei);

        // Ignore error in empty directions
        if ((foundOwn || foundNei) && mag(solutionD & (dr/magdr)) > 0.1)
        {
            scalar dxdr = (x[nei] - x[own])/magdr;
            scalar xc = (x[nei] + x[own])*0.5;
            scalar dl = (dL[own] + dL[nei])*0.5;
            scalar dxDotOwn = gradX[own] & (dr/magdr);
            scalar dxDotNei = gradX[nei] & (-dr/magdr);
            scalar eT =
                Foam::max
                (
                    mag(dxdr - dxDotNei)/(0.3*xc/dl + mag(dxDotNei)),
                    mag(dxdr - dxDotOwn)/(0.3*xc/dl + mag(dxDotOwn))
                );
            if (foundOwn)
            {
                error_[own] = Foam::max(error_[own], eT);
            }
            if (foundNei)
            {
                error_[nei] = Foam::max(error_[nei], eT);
            }
        }
    }

    // Boundary faces
    forAll(error_.boundaryField(), patchi)
    {
        if (error_.boundaryField()[patchi].coupled())
        {
            const fvPatch& patch = x.boundaryField()[patchi].patch();

            const labelUList& faceCells = patch.faceCells();
            scalarField xp
            (
                x.boundaryField()[patchi].patchInternalField()
            );
            scalarField xn
            (
                x.boundaryField()[patchi].patchNeighbourField()
            );
            vectorField drField(patch.delta());
            vectorField gradXp
            (
                gradX.boundaryField()[patchi].patchInternalField()
            );
            vectorField gradXn
            (
                gradX.boundaryField()[patchi].patchNeighbourField()
            );


            forAll(faceCells, facei)
            {
                const vector& dr = drField[facei];
                const scalar magdr = mag(dr);

                // Ignore error in empty directions
                if
                (
                    eCells.found(faceCells[facei])
                 && mag(solutionD & (dr/magdr)) > 0.1
                )
                {
                    scalar dxdr = (xn[facei] - xp[facei])/magdr;
                    scalar xc = (xn[facei] + xp[facei])*0.5;
                    scalar dl = dL[faceCells[facei]];
                    scalar dxDotOwn = gradXp[facei] & (dr/magdr);
                    scalar dxDotNei = gradXn[facei] & (-dr/magdr);
                    scalar eT =
                        Foam::max
                        (
                            mag(dxdr - dxDotNei)
                           /(0.3*xc/dl + mag(dxDotNei)),
                            mag(dxdr - dxDotOwn)
                           /(0.3*xc/dl + mag(dxDotOwn))
                        );
                    error_[faceCells[facei]] =
                        Foam::max(error_[faceCells[facei]], eT);
                }
            }
        }
    }
    if (scale)
    {
        normalize(error_, eCells);
    }
}

// ************************************************************************* //

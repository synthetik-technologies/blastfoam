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
    fieldName_(dict.lookupBackwardsCompatible({"fieldName", "field"}))
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

    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");

    volVectorField gradRho(fvc::grad(rho));
    const scalarField& dL(meshSizeObject::New(mesh_).dx());

    const labelHashSet& eCells = this->errorCells();
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
            scalar dRhodr = (rho[nei] - rho[own])/magdr;
            scalar rhoc = (rho[nei] + rho[own])*0.5;
            scalar dl = (dL[own] + dL[nei])*0.5;
            scalar dRhoDotOwn = gradRho[own] & (dr/magdr);
            scalar dRhoDotNei = gradRho[nei] & (-dr/magdr);
            scalar eT =
                Foam::max
                (
                    mag(dRhodr - dRhoDotNei)/(0.3*rhoc/dl + mag(dRhoDotNei)),
                    mag(dRhodr - dRhoDotOwn)/(0.3*rhoc/dl + mag(dRhoDotOwn))
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
            const fvPatch& patch = rho.boundaryField()[patchi].patch();

            const labelUList& faceCells = patch.faceCells();
            scalarField rhop
            (
                rho.boundaryField()[patchi].patchInternalField()
            );
            scalarField rhon
            (
                rho.boundaryField()[patchi].patchNeighbourField()
            );
            vectorField drField(patch.delta());
            vectorField gradRhop
            (
                gradRho.boundaryField()[patchi].patchInternalField()
            );
            vectorField gradRhon
            (
                gradRho.boundaryField()[patchi].patchNeighbourField()
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
                    scalar dRhodr = (rhon[facei] - rhop[facei])/magdr;
                    scalar rhoc = (rhon[facei] + rhop[facei])*0.5;
                    scalar dl = dL[faceCells[facei]];
                    scalar dRhoDotOwn = gradRhop[facei] & (dr/magdr);
                    scalar dRhoDotNei = gradRhon[facei] & (-dr/magdr);
                    scalar eT =
                        Foam::max
                        (
                            mag(dRhodr - dRhoDotNei)
                           /(0.3*rhoc/dl + mag(dRhoDotNei)),
                            mag(dRhodr - dRhoDotOwn)
                           /(0.3*rhoc/dl + mag(dRhoDotOwn))
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

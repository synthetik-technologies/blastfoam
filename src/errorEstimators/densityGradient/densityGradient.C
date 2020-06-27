/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "densityGradient.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace errorEstimators
{
    defineTypeNameAndDebug(densityGradient, 0);
    addToRunTimeSelectionTable(errorEstimator, densityGradient, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorEstimators::densityGradient::densityGradient
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    errorEstimator(mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimators::densityGradient::~densityGradient()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::errorEstimators::densityGradient::update()
{
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    volScalarField& error(*this);

    volVectorField gradRho(fvc::grad(rho));
    scalarField dL(mesh_.V()/fvc::surfaceSum(mesh_.magSf()));

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const label nInternalFaces = mesh_.nInternalFaces();
    error = 0.0;

    vector solutionD((vector(mesh_.geometricD()) + vector::one)/2.0);

    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        label own = owner[facei];
        label nei = neighbour[facei];
        vector dr = mesh_.C()[nei] - mesh_.C()[own];
        scalar magdr = mag(dr);

        // Ignore error in empty directions
        if (mag(solutionD & (dr/magdr)) > 0.1)
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
            error[own] = Foam::max(error[own], eT);
            error[nei] = Foam::max(error[nei], eT);
        }
    }

    // Boundary faces
    forAll(error.boundaryField(), patchi)
    {
        if (error.boundaryField()[patchi].coupled())
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
                vector dr = drField[facei];
                scalar magdr = mag(dr);

                // Ignore error in empty directions
                if (mag(solutionD & (dr/magdr)) > 0.1)
                {
                    scalar dRhodr = (rhon[facei] - rhop[facei])/magdr;
                    scalar rhoc = (rhon[facei] + rhop[facei])*0.5;
                    scalar dl = dL[faceCells[facei]];
                    scalar dRhoDotOwn = gradRhop[facei] & (dr/magdr);
                    scalar dRhoDotNei = gradRhon[facei] & (-dr/magdr);
                    scalar eT =
                        Foam::max
                        (
                            mag(dRhodr - dRhoDotNei)/(0.3*rhoc/dl + mag(dRhoDotNei)),
                            mag(dRhodr - dRhoDotOwn)/(0.3*rhoc/dl + mag(dRhoDotOwn))
                        );
                    error[faceCells[facei]] =
                        Foam::max(error[faceCells[facei]], eT);
                }
            }
        }
    }
    error.correctBoundaryConditions();
}

// ************************************************************************* //

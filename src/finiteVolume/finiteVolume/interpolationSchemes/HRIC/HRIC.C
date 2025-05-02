/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "HRIC.H"
#include "fvc.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "upwind.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::HRIC::phiface
(
    const scalar cdWeight,
    const scalar faceFlux,
    const scalar phiP,
    const scalar phiN,
    const vector& gradcP,
    const vector& gradcN,
    const vector& d
) const
{
    // Note: 0-1 bounds stabilisation is done in weights calculation
    // HJ, 23/Nov/2011

    // Calculate upwind value, faceFlux C tilde and do stabilisation

    // Upwind, central and downwind value
    scalar phiupw = 0;
    scalar phidown = 0;

    scalar phict = 0;
    scalar cosTheta = 0;

    if (faceFlux > 0)
    {
        cosTheta = mag((gradcP & d)/stabilise(mag(gradcP)*mag(d), small));

        phiupw = phiN - 2.0*(gradcP & d);
        phidown = phiN;

        phiupw = max(min(phiupw, 1), 0);

        phict = (phiP - phiupw)/stabilise(phiN - phiupw, small);
    }
    else
    {
        cosTheta = mag((gradcN & d)/(mag(gradcN)*mag(d) + SMALL));

        phiupw = phiP + 2.0*(gradcN & d);
        phidown = phiP;

        phiupw = max(min(phiupw, 1), 0);

        phict = (phiN - phiupw)/stabilise(phiP - phiupw, small);
    }

    // Calculate phift without Courant number correction
    scalar phift = 0;

    if (phict > 0 && phict <= 0.5)         // use linear upwind scheme
    {
        phift = 2.0*phict;
    }
    else if (phict > 0.5 && phict < 1)     // use downwind
    {
        phift = 1;
    }
    else                                   // use upwind
    {
        phift = phict;
    }

    // Calculate the weighting factors for HRIC
    scalar sqrtCosTheta = sqrt(cosTheta);

    scalar fact = phift*sqrtCosTheta + phict*(1.0 - sqrtCosTheta);

    // Return face value
    return fact*(phidown - phiupw) + phiupw;
}


Foam::scalar Foam::HRIC::limiter
(
    const scalar cdWeight,
    const scalar faceFlux,
    const scalar phiP,
    const scalar phiN,
    const vector& gradcP,
    const vector& gradcN,
    const scalar Cof,
    const vector& d
) const
{
    // Calculate HRIC face value
    scalar phif = phiface
    (
        cdWeight,
        faceFlux,
        phiP,
        phiN,
        gradcP,
        gradcN,
        d
    );
    // Calculate UD and CD face value
    scalar phiU = faceFlux >= 0 ? phiP : phiN;
    scalar phiCD = cdWeight*phiP + (1.0 - cdWeight)*phiN;

    // Calculate the effective limiter for the HRIC interpolation
    scalar CLimiter = max
    (
        min
        (
            (phif - phiU)/stabilise(phiCD - phiU, SMALL),
            2
        ),
        0
    );

    // Courant number correction
    if (Cof < 0.3)
    {
        // Return limiter without Co correction
        return CLimiter;
    }
    else if (Cof >= 0.3 && Cof < 0.7)
    {
        // Blend limiter from no correction at Co = 0.3 to upwind at Co = 0.7
        return CLimiter*(0.7 - Cof)/(0.7 - 0.3);
    }
    else
    {
        return 0;
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::HRIC::limiter
(
    const volScalarField& phi
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tLimiter
    (
        surfaceScalarField::New
        (
            type() + "Limiter(" + phi.name() + ')',
            mesh,
            dimless
        )
    );
    surfaceScalarField& lim = tLimiter.ref();

    volVectorField gradc(fvc::grad(phi));

    surfaceScalarField Cof
    (
        mesh.time().deltaT()
       *upwind<scalar>(mesh, faceFlux_).interpolate
        (
            fvc::surfaceIntegrate(faceFlux_)
        )
    );

    const surfaceScalarField& CDweights = mesh.surfaceInterpolation::weights();

    const UList<label>& owner = mesh.owner();
    const UList<label>& neighbour = mesh.neighbour();

    const vectorField& C = mesh.C();

    scalarField& pLim = lim.primitiveFieldRef();

    forAll(pLim, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        pLim[faceI] = limiter
        (
            CDweights[faceI],
            this->faceFlux_[faceI],
            phi[own],
            phi[nei],
            gradc[own],
            gradc[nei],
            Cof[faceI],
            C[nei] - C[own]
        );
    }

    surfaceScalarField::Boundary& bLim = lim.boundaryFieldRef();

    forAll(bLim, patchi)
    {
        scalarField& pLim = bLim[patchi];

        if (bLim[patchi].coupled())
        {
            const scalarField& pCDweights = CDweights.boundaryField()[patchi];

            const scalarField& pFaceFlux =
                this->faceFlux_.boundaryField()[patchi];

            scalarField pphiP
            (
                phi.boundaryField()[patchi].patchInternalField()
            );

            scalarField pphiN
            (
                phi.boundaryField()[patchi].patchNeighbourField()
            );

            vectorField pGradcP
            (
                gradc.boundaryField()[patchi].patchInternalField()
            );

            vectorField pGradcN
            (
                gradc.boundaryField()[patchi].patchNeighbourField()
            );

            const scalarField& pCof = Cof.boundaryField()[patchi];

            // Build the d-vectors
            // Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
            vectorField pd(bLim[patchi].patch().delta());

            forAll(pLim, faceI)
            {
                pLim[faceI] = limiter
                (
                    pCDweights[faceI],
                    pFaceFlux[faceI],
                    pphiP[faceI],
                    pphiN[faceI],
                    pGradcP[faceI],
                    pGradcN[faceI],
                    pCof[faceI],
                    pd[faceI]
                );
            }
        }
        else
        {
            pLim = 1;
        }
    }

    return tLimiter;
}


namespace Foam
{
defineTypeNameAndDebug(HRIC, 0);

surfaceInterpolationScheme<scalar>::addMeshConstructorToTable<HRIC>
    addHRICMeshConstructorToTable_;

surfaceInterpolationScheme<scalar>::addMeshFluxConstructorToTable<HRIC>
    addHRICMeshFluxConstructorToTable_;

limitedSurfaceInterpolationScheme<scalar>::addMeshConstructorToTable<HRIC>
    addHRICMeshConstructorToLimitedTable_;

limitedSurfaceInterpolationScheme<scalar>::
addMeshFluxConstructorToTable<HRIC>
    addHRICDMeshFluxConstructorToLimitedTable_;
}

// ************************************************************************* //

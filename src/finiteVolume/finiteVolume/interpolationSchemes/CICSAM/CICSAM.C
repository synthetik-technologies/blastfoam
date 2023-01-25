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

#include "CICSAM.H"
#include "fvc.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "upwind.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::CICSAM::lowerBound_ = 0;
const Foam::scalar Foam::CICSAM::upperBound_ = 1;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::CICSAM::limiter
(
    const scalar cdWeight,
    const scalar faceFlux,
    const scalar& phiP,
    const scalar& phiN,
    const vector& gradcP,
    const vector& gradcN,
    const scalar Cof,
    const vector d
) const
{
    // Additional stabilisation for phi out of bounds
    if
    (
        (faceFlux > 0 && (phiP < lowerBound_ || phiN > upperBound_))
     || (faceFlux < 0 && (phiN < lowerBound_ || phiP > upperBound_))
    )
    {
        return 0;
    }

    // Calculate CICSAM weight
    scalar w = weight
    (
        cdWeight,
        faceFlux,
        phiP,
        phiN,
        gradcP,
        gradcN,
        Cof,
        d
    );

    // Calculate CICSAM face value
    scalar phif = w*phiP + (1 - w)*phiN;

    // Calculate UD and CD face value
    scalar phiU = faceFlux >= 0 ? phiP : phiN;
    scalar phiCD = cdWeight*phiP + (1 - cdWeight)*phiN;

    // Calculate the effective limiter for the CICSAM interpolation
    scalar CLimiter = (phif - phiU)/stabilise(phiCD - phiU, small);

    // Limit the limiter between upwind and downwind
    return max(min(CLimiter, 2), 0);
}


Foam::scalar Foam::CICSAM::weight
(
    const scalar cdWeight,
    const scalar faceFlux,
    const scalar& phiP,
    const scalar& phiN,
    const vector& gradcP,
    const vector& gradcN,
    const scalar Cof,
    const vector d
) const
{
    // Additional 0-1 stabilisation.  HJ, 23/Nov/2011
    if
    (
        (faceFlux > 0 && (phiP <= lowerBound_ || phiN >= upperBound_))
    )
    {
        return 1;
    }
    else if
    (
        (faceFlux < 0 && (phiN <= lowerBound_ || phiP >= upperBound_))
    )
    {
        return 0;
    }


    // Calculate upwind value, faceFlux C tilde and do a stabilisation

    scalar phict = 0;
    scalar phiupw = 0;
    scalar costheta = 0;

    if (faceFlux > 0)
    {
        costheta = mag((gradcP & d)/(mag(gradcP)*mag(d) + small));

        phiupw = phiN - 2.0*(gradcP & d);

        phiupw = max(min(phiupw, 1.0), 0.0);

        if ((phiN - phiupw) > 0)
        {
            phict = (phiP - phiupw)/(phiN - phiupw + small);
        }
        else
        {
            phict = (phiP - phiupw)/(phiN - phiupw - small);
        }
    }
    else
    {
        costheta = mag((gradcN & d)/(mag(gradcN)*mag(d) + small));

        phiupw = phiP + 2.0*(gradcN & d);

        phiupw = max(min(phiupw, 1.0), 0.0);

        if ((phiP - phiupw) > 0)
        {
            phict = (phiN - phiupw)/(phiP - phiupw + small);
        }
        else
        {
            phict = (phiN - phiupw)/(phiP - phiupw - small);
        }
    }


    // Calculate the weighting factors for CICSAM

    scalar cicsamFactor = (k_ + small)/(1 - k_ + small);

    costheta = min(1.0, cicsamFactor*(costheta));
    costheta = (cos(2*(acos(costheta))) + 1.0)/2.0;

    scalar k1 = 3.0*Cof*(Cof - 1.0)/(2.0*Cof*Cof + 6.0*Cof - 8.0);
    scalar k2 = Cof;
    scalar k3 = (3.0*Cof + 5.0)/(2.0*Cof + 6.0);
    scalar weight;

    if (phict == 1)
    {
        weight = phict >= k3 ? 1.0 : 0.0;
    }
    else if (phict > 0 && phict <= k1)             // use blended scheme 1
    {
        scalar phifCM = phict/(Cof + small);
        weight = (phifCM - phict)/(1.0 - phict);
    }
    else if (phict > k1 && phict <= k2)     // use blended scheme 2
    {
        scalar phifHC = phict/(Cof + small);
        scalar phifUQ = (8.0*Cof*phict + (1.0 - Cof)*(6.0*phict + 3.0))/8.0;
        scalar phifCM = costheta*phifHC + (1.0 - costheta)*phifUQ;
        weight = (phifCM - phict)/(1.0 - phict);
    }
    else if (phict > k2 && phict < k3)     // use blended scheme 3
    {
        scalar phifUQ = (8.0*Cof*phict + (1.0 - Cof)*(6.0*phict + 3.0))/8.0;
        scalar phifCM = costheta + (1.0 - costheta)*phifUQ;
        weight = (phifCM - phict)/(1.0 - phict);
    }
    else if (phict >= k3 && phict <= 1)     // use downwind
    {
        weight = 1;
    }
    else                                    // use upwind
    {
        weight = 0;
    }

    // Correction: In original code, Onno's weights were calculated
    // the wrong way around
    // HJ, 24/Nov/2011
    if (faceFlux > 0)
    {
        return 1 - weight;
    }
    else
    {
        return weight;
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::CICSAM::limiter
(
    const volScalarField& phi
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tLimiter
    (
        new surfaceScalarField
        (
            IOobject
            (
                type() + "Limiter(" + phi.name() + ')',
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimless
        )
    );

    //Note: Changing this line may mess up conversion to old API style
    surfaceScalarField& lim = tLimiter.ref();

    volVectorField gradc(fvc::grad(phi));

/*
    surfaceScalarField Cof =
        mesh.time().deltaT()
       *upwind<scalar>(mesh, faceFlux_).interpolate
        (
            fvc::surfaceIntegrate(faceFlux_)
        );
*/

    surfaceScalarField Cof
    (
        0.5 * mesh.time().deltaT()
        *upwind<scalar>(mesh, faceFlux_).interpolate
        (
            fvc::surfaceIntegrate(mag(faceFlux_))
        )
    );


    const surfaceScalarField& CDweights = mesh.surfaceInterpolation::weights();

    const UList<label>& owner = mesh.owner();
    const UList<label>& neighbour = mesh.neighbour();

    const vectorField& C = mesh.C();

    scalarField& pLim = lim.ref();

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
            pLim = 1.0;
        }
    }

    return tLimiter;
}



Foam::tmp<Foam::surfaceScalarField> Foam::CICSAM::weights
(
    const volScalarField& phi
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tWeightingFactors
    (
        new surfaceScalarField(mesh.surfaceInterpolation::weights())
    );
    //Note: Changing this line may mess up conversion to old API style
    surfaceScalarField& weightingFactors = tWeightingFactors.ref();

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

    scalarField& w = weightingFactors.ref();

    forAll(w, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        w[faceI] = weight
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

    surfaceScalarField::Boundary& bWeights =
        weightingFactors.boundaryFieldRef();

    forAll(bWeights, patchi)
    {
        scalarField& pWeights = bWeights[patchi];

        if (bWeights[patchi].coupled())
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
            vectorField pd(bWeights[patchi].patch().delta());

            forAll(pWeights, faceI)
            {
                pWeights[faceI] = weight
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
            pWeights = 1.0;
        }
    }

    return tWeightingFactors;
}


namespace Foam
{
defineTypeNameAndDebug(CICSAM, 0);

surfaceInterpolationScheme<scalar>::addMeshConstructorToTable<CICSAM>
    addCICSAMMeshConstructorToTable_;

surfaceInterpolationScheme<scalar>::addMeshFluxConstructorToTable<CICSAM>
    addCICSAMMeshFluxConstructorToTable_;

limitedSurfaceInterpolationScheme<scalar>::addMeshConstructorToTable<CICSAM>
    addCICSAMMeshConstructorToLimitedTable_;

limitedSurfaceInterpolationScheme<scalar>::
addMeshFluxConstructorToTable<CICSAM>
    addCICSAMMeshFluxConstructorToLimitedTable_;
}

// ************************************************************************* //

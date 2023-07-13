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

#include "mHRIC.H"
#include "fvc.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "upwind.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::mHRIC::phiface
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
    // Upwind cell
    scalar phiU = 0;

    // Donor cell
    scalar phiD = 0;

    // Acceptor cell
    scalar phiA = 0;

    // phiP: face owner phi
    // phiN: face neighbor phi
    // d: distance vector -> rN - rP

    // Nomalised Donor cell phi
    scalar phiDTilda = 0;
    scalar costheta = 0;


    //        ?----P--f--N----   :u>0
    //        ----U----D--f--A----

    //         ----P--f--N----?   :u<0
    //        ---------A--f--D----U
    if (faceFlux > 0)
    {
        phiD = min(max(phiP, 0), 1);
        phiA = min(max(phiN, 0), 1);
        phiU = phiA - 2.0*(gradcP & d);
        phiU = min(max(phiU, 0), 1);

        if (phiA - phiU < small)
        {
            phiDTilda = 0;
        }
        else
        {
            if (phiA - phiU >= 0)   // phiU < phiD < phiA
            {
                phiD = max(min(phiD, phiA), phiU);
            }
            else                    // phiU > phiD > phiA
            {
                phiD = max(min(phiD, phiU), phiA);
            }
            phiDTilda = (phiD - phiU)/(phiA - phiU);
        }

        if (mag(gradcP) < small)
        {
            costheta = 0;
        }
        else
        {
            costheta = (gradcP & d)/(mag(gradcP)*mag(d));
        }
    }
    else
    {
        phiD = min(max(phiN, 0), 1);
        phiA = min(max(phiP, 0), 1);
        phiU = phiA + 2.0*(gradcN & d);
        phiU = min(max(phiU, 0), 1);

        if (phiA - phiU < small)
        {
            phiDTilda = 0;
        }
        else
        {
            if (phiA - phiU >= 0)   // phiU < phiD < phiA
            {
                phiD = max(min(phiD, phiA), phiU);
            }
            else                    // phiU > phiD > phiA
            {
                phiD = max(min(phiD, phiU), phiA);
            }
            phiDTilda = (phiD - phiU)/(phiA - phiU);
        }
        if (mag(gradcN) < small)
        {
            costheta = 0;
        }
        else
        {
            costheta = (gradcN & d)/(mag(gradcN)*mag(d));
        }
    }

    phiDTilda = min(max(mag(phiDTilda), 0), 1);

    // phiFTilda
    scalar phiFTilda = 0;

    if (phiDTilda >= 0 && phiDTilda <= 0.5)
    {
        phiFTilda = 2.0*phiDTilda;
    }
    else if (phiDTilda >= 0.5 && phiDTilda <= 1)
    {
        phiFTilda = 1;
    }
    else
    {
        phiFTilda = phiDTilda;
    }

    // phifUQTilda
    scalar phifUQTilda = 0;

    if (phiDTilda >= 0 && phiDTilda <= 1)
    {
        phifUQTilda = min(phiFTilda, (6.0*phiDTilda + 3.0)/8.0);
    }
    else
    {
        phifUQTilda = phiDTilda;
    }

    scalar gammarf = sqrt(mag(costheta));

    // phiFsTilda
    scalar phiFsTilda = gammarf*phiFTilda + (1.0 - gammarf)*phifUQTilda;

    // Return face value
    return phiFsTilda*(phiA - phiU) + phiA;
}


Foam::scalar Foam::mHRIC::limiter
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
    // Calculate mHRIC face value
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

    // Calculate the effective limiter for the mHRIC interpolation
    scalar CLimiter = max
    (
        min
        (
            (phif - phiU)/stabilise(phiCD - phiU, small),
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


Foam::tmp<Foam::surfaceScalarField> Foam::mHRIC::limiter
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
            pLim = 1;
        }
    }

    return tLimiter;
}


namespace Foam
{
defineTypeNameAndDebug(mHRIC, 0);

surfaceInterpolationScheme<scalar>::addMeshConstructorToTable<mHRIC>
    addmHRICMeshConstructorToTable_;

surfaceInterpolationScheme<scalar>::addMeshFluxConstructorToTable<mHRIC>
    addmHRICMeshFluxConstructorToTable_;

limitedSurfaceInterpolationScheme<scalar>::addMeshConstructorToTable<mHRIC>
    addmHRICMeshConstructorToLimitedTable_;

limitedSurfaceInterpolationScheme<scalar>::
addMeshFluxConstructorToTable<mHRIC>
    addmHRICDMeshFluxConstructorToLimitedTable_;
}

// ************************************************************************* //

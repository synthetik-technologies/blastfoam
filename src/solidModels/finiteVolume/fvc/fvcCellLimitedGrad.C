/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvcCellLimitedGrad.H"

#include "Field.H"
#include "primitivePatch.H"
#include "fvc.H"

#include "fvcGradf.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"

#include "gaussGrad.H"
#ifdef FOAMEXTEND
    #include "extendedLeastSquaresGrad.H"
#endif
#include "leastSquaresGrad.H"
//#include "leastSquaresBiLinearGrad.H"

#include "volMesh.H"
#include "surfaceMesh.H"
#include "fixedValueFvPatchFields.H"
#include "IStringStream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


inline void limitFace
(
    scalar& limiter,
    const scalar& maxDelta,
    const scalar& minDelta,
    const scalar& extrapolate
)
{
    if (extrapolate > maxDelta + VSMALL)
    {
        limiter = min(limiter, maxDelta/extrapolate);
    }
    else if (extrapolate < minDelta - VSMALL)
    {
        limiter = min(limiter, minDelta/extrapolate);
    }
}

inline void limitFace
(
    vector& limiter,
    const vector& maxDelta,
    const vector& minDelta,
    const vector& extrapolate
)
{
    for(direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        limitFace
        (
            limiter.component(cmpt),
            maxDelta.component(cmpt),
            minDelta.component(cmpt),
            extrapolate.component(cmpt)
        );
    }
}


tmp<volVectorField> cellLimitedGrad
(
    const volScalarField& vsf,
    const pointScalarField& pf,
    const bool lsCellGrad,
    const scalar k
)
{
    const fvMesh& mesh = vsf.mesh();

    tmp<GeometricField<vector, fvPatchField, volMesh> > tGrad
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "grad(" + vsf.name() + ")",
                vsf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<vector>
            (
                "0",
                vsf.dimensions()/dimLength,
                pTraits<vector>::zero
            ),
            zeroGradientFvPatchField<vector>::typeName
        )
    );

#ifdef FOAMEXTEND
    if (lsCellGrad)
    {
        // Extended least squares
        IStringStream stream("0");
        tGrad =
            fv::extendedLeastSquaresGrad<scalar>
            (
                mesh,
                stream
            ).grad(vsf);

        // Gauss grad
//         tGrad = fv::gaussGrad<scalar>(mesh).grad(vsf);
    }
    else
    {
        tGrad = fvc::grad(vsf, pf);
    }

//     tmp<volVectorField> tGrad = basicGradScheme_().grad(vsf);

    if (k < SMALL)
    {
        return tGrad;
    }

    volVectorField& g = tGrad();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    scalarField maxVsf(vsf.internalField());
    scalarField minVsf(vsf.internalField());

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        scalar vsfOwn = vsf[own];
        scalar vsfNei = vsf[nei];

        maxVsf[own] = max(maxVsf[own], vsfNei);
        minVsf[own] = min(minVsf[own], vsfNei);

        maxVsf[nei] = max(maxVsf[nei], vsfOwn);
        minVsf[nei] = min(minVsf[nei], vsfOwn);
    }


    const volScalarField::GeometricBoundaryField& bsf = vsf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvPatchScalarField& psf = bsf[patchi];

        const unallocLabelList& pOwner = mesh.boundary()[patchi].faceCells();

        if (psf.coupled())
        {
            scalarField psfNei = psf.patchNeighbourField();

            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                scalar vsfNei = psfNei[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
        else
        {
            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                scalar vsfNei = psf[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
    }

    maxVsf -= vsf;
    minVsf -= vsf;

    if (k < 1.0)
    {
        scalarField maxMinVsf = (1.0/k - 1.0)*(maxVsf - minVsf);
        maxVsf += maxMinVsf;
        minVsf -= maxMinVsf;

        //maxVsf *= 1.0/k_;
        //minVsf *= 1.0/k_;
    }


    // create limiter
    scalarField limiter(vsf.internalField().size(), 1.0);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // owner side
        limitFace
        (
            limiter[own],
            maxVsf[own],
            minVsf[own],
            (Cf[facei] - C[own]) & g[own]
        );

        // neighbour side
        limitFace
        (
            limiter[nei],
            maxVsf[nei],
            minVsf[nei],
            (Cf[facei] - C[nei]) & g[nei]
        );
    }

    forAll(bsf, patchi)
    {
        const unallocLabelList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            label own = pOwner[pFacei];

            limitFace
            (
                limiter[own],
                maxVsf[own],
                minVsf[own],
                (pCf[pFacei] - C[own]) & g[own]
            );
        }
    }

    if (fv::debug)
    {
        Info<< "gradient limiter for: " << vsf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    g.internalField() *= limiter;
    g.correctBoundaryConditions();
    fv::gaussGrad<scalar>::correctBoundaryConditions(vsf, g);
#else
    notImplemented("Not implemented for this version of OpenFOAM/FOAM");
#endif

    return tGrad;
}


tmp<volTensorField> cellLimitedGrad
(
    const volVectorField& vf,
    const pointVectorField& pf,
    const bool lsCellGrad,
    const scalar k
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<GeometricField<tensor, fvPatchField, volMesh> > tGrad
    (
        new GeometricField<tensor, fvPatchField, volMesh>
        (
            IOobject
            (
                "grad(" + vf.name() + ")",
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<tensor>
            (
                "0",
                vf.dimensions()/dimLength,
                pTraits<tensor>::zero
            ),
            zeroGradientFvPatchField<tensor>::typeName
        )
    );

#ifdef FOAMEXTEND
    if (lsCellGrad)
    {
        // // Least squares
        // tGrad = fv::leastSquaresBiLinearGrad<vector>(mesh).grad(vf);

        IStringStream stream("0");
        tGrad =
            fv::extendedLeastSquaresGrad<vector>
            (
                mesh,
                stream
            ).grad(vf);
    }
    else
    {
        tGrad = fvc::grad(vf, pf);
    }

//     tmp<volTensorField> tGrad = fvc::grad(vf, pf);

    if (k < SMALL)
    {
        return tGrad;
    }

    volTensorField& g = tGrad();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    vectorField maxVsf(vf.internalField());
    vectorField minVsf(vf.internalField());

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        const vector& vsfOwn = vf[own];
        const vector& vsfNei = vf[nei];

        maxVsf[own] = max(maxVsf[own], vsfNei);
        minVsf[own] = min(minVsf[own], vsfNei);

        maxVsf[nei] = max(maxVsf[nei], vsfOwn);
        minVsf[nei] = min(minVsf[nei], vsfOwn);
    }

    const volVectorField::GeometricBoundaryField& bsf = vf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvPatchVectorField& psf = bsf[patchi];
        const unallocLabelList& pOwner = mesh.boundary()[patchi].faceCells();

        if (psf.coupled())
        {
            vectorField psfNei = psf.patchNeighbourField();

            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                const vector& vsfNei = psfNei[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
        else
        {
            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                const vector& vsfNei = psf[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
    }

    maxVsf -= vf;
    minVsf -= vf;

    if (k < 1.0)
    {
        vectorField maxMinVsf = (1.0/k - 1.0)*(maxVsf - minVsf);
        maxVsf += maxMinVsf;
        minVsf -= maxMinVsf;

        //maxVsf *= 1.0/k_;
        //minVsf *= 1.0/k_;
    }


    // create limiter
    vectorField limiter(vf.internalField().size(), vector::one);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // owner side
        limitFace
        (
            limiter[own],
            maxVsf[own],
            minVsf[own],
            (Cf[facei] - C[own]) & g[own]
        );

        // neighbour side
        limitFace
        (
            limiter[nei],
            maxVsf[nei],
            minVsf[nei],
            (Cf[facei] - C[nei]) & g[nei]
        );
    }

    forAll(bsf, patchi)
    {
        const unallocLabelList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            label own = pOwner[pFacei];

            limitFace
            (
                limiter[own],
                maxVsf[own],
                minVsf[own],
                ((pCf[pFacei] - C[own]) & g[own])
            );
        }
    }

    if (false)
    {
        Info<< "gradient limiter for: " << vf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    tensorField& gIf = g.internalField();

    forAll(gIf, celli)
    {
        gIf[celli] = tensor
        (
            cmptMultiply(limiter[celli], gIf[celli].x()),
            cmptMultiply(limiter[celli], gIf[celli].y()),
            cmptMultiply(limiter[celli], gIf[celli].z())
        );
    }

    g.correctBoundaryConditions();
    fv::gaussGrad<vector>::correctBoundaryConditions(vf, g);
#else
    notImplemented("Not implemented for this version of OpenFOAM/FOAM");
#endif

    return tGrad;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

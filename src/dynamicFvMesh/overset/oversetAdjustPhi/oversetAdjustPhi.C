/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "oversetAdjustPhi.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cellCellStencilObject.H"
#include "syncTools.H"
#include "fv.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::oversetAdjustPhi
(
    surfaceScalarField& phi,
    const volVectorField& U
)
{
    const fvMesh& mesh = U.mesh();

    const cellCellStencilObject& overlap = Stencil::New(mesh);
    const labelList& cellTypes = overlap.cellTypes();
    const labelList& zoneID = overlap.zoneID();
    label nZones = gMax(zoneID)+1;



    // Pass1: accumulate all fluxes, calculate correction factor

    scalarField massIn(nZones, Zero);
    scalarField adjustableMassOut(nZones, Zero);

    surfaceScalarField::Boundary& bphi =
        phi.boundaryFieldRef();


    // Check all faces on the outside of interpolated cells
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    {
        forAll(own, facei)
        {
            label zonei = zoneID[own[facei]];   // note:own and nei in same zone

            label ownType = cellTypes[own[facei]];
            label neiType = cellTypes[nei[facei]];

            bool ownCalc =
                (ownType == cellCellStencil::CALCULATED)
             && (neiType == cellCellStencil::INTERPOLATED);

            bool neiCalc =
                (ownType == cellCellStencil::INTERPOLATED)
             && (neiType == cellCellStencil::CALCULATED);


            if (ownCalc || neiCalc)
            {
                // Calculate flux w.r.t. calculated cell
                scalar flux = phi[facei];
                if (neiCalc)
                {
                    flux = -flux;
                }

                if (flux < 0.0)
                {
                    massIn[zonei] -= flux;
                }
                else
                {
                    adjustableMassOut[zonei] += flux;
                }
            }
        }
    }


    // Check all coupled faces on the outside of interpolated cells
    labelList neiCellTypes;
    syncTools::swapBoundaryCellList(mesh, cellTypes, neiCellTypes);
    {
        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            const fvsPatchScalarField& phip = bphi[patchi];
            const labelUList& fc = Up.patch().faceCells();

            label start = Up.patch().start();

            forAll(fc, i)
            {
                label facei = start+i;
                label celli = fc[i];
                label ownType = cellTypes[celli];
                label neiType = neiCellTypes[facei-mesh.nInternalFaces()];

                bool ownCalc =
                    (ownType == cellCellStencil::CALCULATED)
                 && (neiType == cellCellStencil::INTERPOLATED);


                if (ownCalc)
                {
                    // Calculate flux w.r.t. calculated cell
                    scalar flux = phip[i];

                    if (flux < 0.0)
                    {
                        massIn[zoneID[celli]] -= flux;
                    }
                    else
                    {
                        adjustableMassOut[zoneID[celli]] += flux;
                    }
                }
            }
        }
    }

    // Calculate the total flux in the domain, used for normalisation
    scalar totalFlux = VSMALL + sum(mag(phi)).value();

    forAll(massIn, zonei)
    {
        reduce(massIn[zonei], sumOp<scalar>());
        reduce(adjustableMassOut[zonei], sumOp<scalar>());
    }


    scalarField massCorr(nZones, 1.0);

    forAll(massIn, zonei)
    {
        scalar magAdjustableMassOut = mag(adjustableMassOut[zonei]);

        if
        (
            magAdjustableMassOut > VSMALL
         && magAdjustableMassOut/totalFlux > SMALL
        )
        {
            massCorr[zonei] = massIn[zonei]/adjustableMassOut[zonei];
        }
        else if (mag(massIn[zonei])/totalFlux > 1e-8)
        {
            WarningInFunction
                << "Continuity error cannot be removed by adjusting the"
                   " flow at fringe faces.\n    Please check the cell types"
                << " from the overset analysis."
                << nl
                << "Zone                    : " << zonei << nl
                << "Total flux              : " << totalFlux << nl
                << "Specified mass inflow   : " << massIn[zonei] << nl
                << "Adjustable mass outflow : " << adjustableMassOut[zonei]
                << nl << endl;
        }


        if (fv::debug)
        {
            Info<< "Zone                    : " << zonei << nl
                << "Total flux              : " << totalFlux << nl
                << "Specified mass inflow   : " << massIn[zonei] << nl
                << "Adjustable mass outflow : " << adjustableMassOut[zonei]
                << nl
                << "Correction factor       : " << massCorr[zonei] << nl
                << endl;
        }
    }



    // Pass2: adjust fluxes

    forAll(own, facei)
    {
        label zonei = zoneID[own[facei]];   // note:own and nei in same zone

        label ownType = cellTypes[own[facei]];
        label neiType = cellTypes[nei[facei]];

        bool ownCalc =
            (ownType == cellCellStencil::CALCULATED)
         && (neiType == cellCellStencil::INTERPOLATED);

        bool neiCalc =
            (ownType == cellCellStencil::INTERPOLATED)
         && (neiType == cellCellStencil::CALCULATED);

        if (ownCalc || neiCalc)
        {
            // Calculate flux w.r.t. calculated cell
            scalar flux = phi[facei];
            if (neiCalc)
            {
                flux = -flux;
            }

            if (flux < 0.0)
            {
                phi[facei] /= Foam::sqrt(massCorr[zonei]);
            }
            else
            {
                phi[facei] *= Foam::sqrt(massCorr[zonei]);
            }
        }
    }

    forAll(bphi, patchi)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        fvsPatchScalarField& phip = bphi[patchi];
        const labelUList& fc = Up.patch().faceCells();

        label start = Up.patch().start();

        forAll(fc, i)
        {
            label facei = start+i;
            label celli = fc[i];
            label zonei = zoneID[celli];   // note:own and nei in same zone
            label ownType = cellTypes[celli];
            label neiType = neiCellTypes[facei-mesh.nInternalFaces()];

            bool ownCalc =
                (ownType == cellCellStencil::CALCULATED)
             && (neiType == cellCellStencil::INTERPOLATED);

            bool neiCalc =
                (ownType == cellCellStencil::INTERPOLATED)
             && (neiType == cellCellStencil::CALCULATED);

            if (ownCalc || neiCalc)
            {
                // Calculate flux w.r.t. calculated cell
                scalar flux = phip[i];
                if (neiCalc)
                {
                    flux = -flux;
                }

                if (flux < 0.0)
                {
                    phip[i] /= Foam::sqrt(massCorr[zonei]);
                }
                else
                {
                    phip[i] *= Foam::sqrt(massCorr[zonei]);
                }
            }
        }
    }

    return sum(mag(massIn))/totalFlux < SMALL
        && sum(mag(adjustableMassOut))/totalFlux < SMALL;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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

#include "findRefCell.H"
#include "regionSplit.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::setRefCells
(
    const volScalarField& field,
    const volScalarField& fieldRef,
    const dictionary& dict,
    boolList& regionNeedReference,
    labelList& regionRefCells,
    scalarField& regionRefValues,
    const bool forceReference
)
{
    const regionSplit& regions = regionSplit::New(field.mesh());

    regionNeedReference.setSize(regions.nRegions(), true);
    regionRefCells.setSize(regions.nRegions(), -1);
    regionRefValues.setSize(regions.nRegions(), 0);

    if (!forceReference)
    {
        const volScalarField::GeometricBoundaryField& bfld =
            fieldRef.boundaryField();

        forAll(bfld, patchI)
        {
            if (bfld[patchI].fixesValue())
            {
                // Unmark all regions

                const labelUList& fc = bfld[patchI].patch().patch().faceCells();

                forAll(fc, faceI)
                {
                    regionNeedReference[regions[fc[faceI]]] = false;
                }
            }
        }

        Pstream::listCombineGather(regionNeedReference, orEqOp<bool>());
        Pstream::listCombineScatter(regionNeedReference);
    }


    label nRefs = 0;
    forAll(regionNeedReference, regionI)
    {
        if (regionNeedReference[regionI])
        {
            nRefs++;
        }
    }

    if (nRefs == 0)
    {
        return;
    }

    // Get the reference cells for all the regions

    word refCellName = field.name() + "RefCells";
    word refPointName = field.name() + "RefPoints";
    word refValueName = field.name() + "RefValues";


    // (per region!) does region have reference cell?
    labelList hasRef(regionNeedReference.size(), Zero);


    const labelList refValues(dict.lookup(refValueName));


    if (dict.found(refCellName))
    {
        // Have specified reference cells (on master!)

        if (Pstream::master())
        {
            labelList refCells(dict.lookup(refCellName));

            if (refCells.size() != regionNeedReference.size())
            {
                FatalIOErrorInFunction(dict)
                    << "Number of refCells " << refCells.size()
                    << " does not correspond to number of regions "
                    << regionNeedReference.size()
                    << exit(FatalIOError);
            }

            forAll(refCells, i)
            {
                label regionI = regions[refCells[i]];

                if (regionNeedReference[regionI])
                {
                    regionRefCells[regionI] = refCells[i];
                    regionRefValues[regionI] = refValues[i];
                }
            }


            forAll(regionNeedReference, regionI)
            {
                if
                (
                    regionNeedReference[regionI]
                 && regionRefCells[regionI] == -1
                )
                {
                    FatalIOErrorInFunction(dict)
                        << "Have no reference cell for region " << regionI
                        << nl
                        << "Overall per-region reference cells "
                        << regionRefCells
                        << exit(FatalIOError);
                }
            }
        }
    }
    else if (dict.found(refPointName))
    {
        pointField refPoints(dict.lookup(refPointName));

        if (refPoints.size() != regionNeedReference.size())
        {
            FatalIOErrorInFunction(dict)
                << "Number of refPoints " << refPoints.size()
                << " does not correspond to number of regions "
                << regionNeedReference.size()
                << exit(FatalIOError);
        }

        labelList hasRef(refPoints.size(), Zero);

        forAll(refPoints, i)
        {
            // Note: find reference cell using facePlanes to avoid constructing
            //       face decomposition structure. Most likely the reference
            //       cell is an undistorted one so this should not be a
            //       problem.

            label celli = field.mesh().findCell
            (
                refPoints[i],
                polyMesh::FACEPLANES
            );

            if (celli >= 0)
            {
                Pout<< "Found point " << refPoints[i]
                    << " in reference cell " << celli
                    << " at " << field.mesh().cellCentres()[celli]
                    << " for region " << regions[celli]
                    << endl;

                regionRefCells[regions[celli]] = celli;
                hasRef[regions[celli]] = 1;
            }
        }

        Pstream::listCombineGather(hasRef, plusEqOp<label>());
        Pstream::listCombineScatter(hasRef);

        forAll(hasRef, regionI)
        {
            if (hasRef[regionI] != 1)
            {
                FatalIOErrorInFunction(dict)
                    << "Unable to set reference cell for field " << field.name()
                    << nl
                    << "    Reference points " << refPointName
                    << " " << refPoints
                    << nl << "    For region " << regionI
                    << " found on " << hasRef[regionI]
                    << " domains (should be one)"
                    << nl << exit(FatalIOError);
            }
        }
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Unable to set reference cell for field " << field.name()
            << nl
            << "    Please supply either " << refCellName
            << " or " << refPointName << nl << exit(FatalIOError);
    }
}


// ************************************************************************* //

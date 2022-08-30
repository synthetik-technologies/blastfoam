/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "setCellDisplacements.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setCellDisplacements, 0);
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::setCellDisplacements::readDict()
{
    Info<< type() << ": reading cellDisplacements" << endl;

    // Lookup the names of each cell displacement
    const wordList cellDispNames = dict_.toc();

    // Initialise fields
    cellIDs_.setSize(cellDispNames.size(), -1);
    constantDisps_.setSize(cellDispNames.size(), vector::zero);
    timeVaryingDisps_.setSize(cellDispNames.size());
    currentCellDisps_.setSize(cellDispNames.size(), vector::zero);

    // Initialise settings for each cell
    forAll(cellDispNames, cellNameI)
    {
        const dictionary& curCellDict =
            dict_.subDict(cellDispNames[cellNameI]);

        // Check if displacements are time-varying or constant
        if (curCellDict.found("timeVaryingDisplacement"))
        {
            timeVaryingDisps_.set
            (
                cellNameI,
                new Function1s::Table<vector>
                (
                    "timeVaryingDisplacement",
                    curCellDict
                )
            );
        }
        else
        {
            constantDisps_[cellNameI] = vector
            (
                curCellDict.lookup("displacement")
            );
        }

        // Lookup the approximate cell coordinate
        const vector coord = vector
        (
            curCellDict.lookup("approximateCoordinate")
        );

        // Find the closest cells in the mesh
        // This cell should only exist on one processor

        scalar dist = GREAT;
        const vectorField& C = mesh_.C().internalField();

        forAll(C, cellI)
        {
            const scalar newDist = mag(C[cellI] - coord);

            if (newDist < dist)
            {
                dist = newDist;
                cellIDs_[cellNameI] = cellI;
            }
        }

        // Find the closest cell globally
        const scalar minDist = returnReduce(dist, minOp<scalar>());
        label procID = int(GREAT);
        if (mag(minDist - dist) > SMALL)
        {
            // -1 signifies that the current proc does not have the closest
            // cell
            cellIDs_[cellNameI] = -1;
        }
        else
        {
            procID = Pstream::myProcNo();
        }

        // If there is more than one processor with the closest cell (could
        // happen if point is specified on a process patch) then we will take
        // the processor with the lowest proc number
        if (Pstream::myProcNo() != returnReduce(procID, minOp<int>()))
        {
            cellIDs_[cellNameI] = -1;
        }

        if (cellIDs_[cellNameI] != -1)
        {
            Pout<< type() << ": desired coordinate = " << coord
                << ", using cell " <<  cellIDs_[cellNameI]
                << " with cell-centre = " << C[cellIDs_[cellNameI]] << endl;
        }
    }

    // In parallel, we will remove cellIDs not on this proc
    if (Pstream::parRun())
    {
        // Count cellIDs that are on this proc (index is not -1)
        int i = 0;
        forAll(cellIDs_, cI)
        {
            if (cellIDs_[cI] != -1)
            {
                i++;
            }
        }

        // Create new lists
        labelList newCellIDs(i, label(-1));
        vectorField newConstantDisps(i, vector::zero);
        PtrList< Function1s::Table<vector> > newTimeVaryingDisps(i);
        vectorField newCurrentCellDisps(i, vector::zero);

        // Copy over cell data
        i = 0;
        forAll(cellIDs_, cI)
        {
            if (cellIDs_[cI] != -1)
            {
                newCellIDs[i] = cellIDs_[cI];
                newConstantDisps[i] = constantDisps_[cI];
                newTimeVaryingDisps.set
                (
                    i,
                    new Function1s::Table<vector>(timeVaryingDisps_[cI])
                );
                newCurrentCellDisps[i] = currentCellDisps_[cI];
                i++;
            }
        }

        // Reset fields
        cellIDs_ = newCellIDs;
        constantDisps_ = newConstantDisps;
        timeVaryingDisps_.transfer(newTimeVaryingDisps);
        currentCellDisps_ = newCurrentCellDisps;

        Pout<< type() << ": proc " << Pstream::myProcNo() << " has "
            << cellIDs_.size() << " cells with setDisplacements" << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setCellDisplacements::setCellDisplacements
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    cellIDs_(),
    constantDisps_(),
    timeVaryingDisps_(),
    currentCellDisps_(),
    curTimeIndex_(-1)
{
    readDict();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField& Foam::setCellDisplacements::cellDisps() const
{
    // Check if it is a new time-step
    if (curTimeIndex_ != mesh_.time().timeIndex())
    {
        curTimeIndex_ = mesh_.time().timeIndex();

        // Update cell displacements
        forAll(cellIDs_, cI)
        {
            const label curCellID = cellIDs_[cI];

            if (curCellID != -1)
            {
                if (timeVaryingDisps_.set(cI))
                {
                    // Time-varying
                    currentCellDisps_[cI] = timeVaryingDisps_[cI].value
                    (
                        mesh_.time().timeOutputValue()
                    );
                }
                else
                {
                    // Constant in time
                    currentCellDisps_[cI] = constantDisps_[cI];
                }
            }
        }
    }

    return currentCellDisps_;
}



// ************************************************************************* //

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
    PtrList<entry> cellDisps(dict_.lookup("cellDisplacements"));

    // Initialise fields
    cellIDs_.setSize(cellDisps.size(), -1);
    disps_.setSize(cellDisps.size());

    // Initialise settings for each cell
    forAll(cellDisps, cellDispI)
    {
        const dictionary& curCellDict = cellDisps[cellDispI].dict();

        // Check if displacements are time-varying or constant
        disps_.set
        (
            cellDispI,
            Function1<vector>::New
            (
                "displacement",
                curCellDict
            ).ptr()
        );

        // Lookup the approximate cell coordinate
        vector coord(curCellDict.lookup("point"));

        // Find the closest cells in the mesh
        // This cell should only exist on one processor
        label cellI = mesh_.findNearestCell(coord);
        scalar error = mag(mesh_.C()[cellI] - coord);
        scalar minError = returnReduce(error, minOp<scalar>());
        if (error != minError)
        {
            cellI = -1;
            coord = -great*vector::one;
        }
        reduce(coord, maxOp<vector>());

        if (cellI != -1)
        {
            Pout<< type() << ": desired coordinate = " << coord
                << ", using cell " <<  cellI
                << " with cell-centre = " << mesh_.C()[cellI] << endl;
        }
        cellIDs_[cellDispI] = cellI;
    }

    // Initialise settings for each cell
    label cellDispI = 0;
    forAll(cellDisps, cellDispi)
    {
        const label cellI = cellIDs_[cellDispi];
        if (cellI != -1)
        {
            if (cellDispI != cellDispi)
            {
                disps_.set(cellDispI, disps_[cellDispi].clone().ptr());
                cellIDs_[cellDispI] = cellI;
            }
            cellDispI++;
        }
    }
    cellIDs_.setSize(cellDispI);
    disps_.setSize(cellDispI);
    currentCellDisps_.setSize(cellDispI, vector::zero);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setCellDisplacements::setCellDisplacements(const fvMesh& mesh)
:
    mesh_(mesh),
    dict_(),
    cellIDs_(),
    disps_(),
    currentCellDisps_(),
    curTimeIndex_(-1)
{}

Foam::setCellDisplacements::setCellDisplacements
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    cellIDs_(),
    disps_(),
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
                // Time-varying
                currentCellDisps_[cI] = disps_[cI].value
                (
                    mesh_.time().timeOutputValue()
                );
            }
        }
    }

    return currentCellDisps_;
}



// ************************************************************************* //

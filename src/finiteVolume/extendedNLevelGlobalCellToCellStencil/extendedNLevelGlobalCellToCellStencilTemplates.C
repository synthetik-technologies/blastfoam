/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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

#include "extendedNLevelGlobalCellToCellStencil.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
void Foam::extendedNLevelGlobalCellToCellStencil::collectData
(
    const Field<Type>& fld,
    List<List<Type>>& stencilFld
) const
{
    // 1. Construct cell data in compact addressing
    List<Type> flatFld(mapPtr_->constructSize(), Zero);

    // Insert my internal values
    forAll(fld, celli)
    {
        flatFld[celli] = fld[celli];
    }

    // Do all swapping
    mapPtr_->distribute(flatFld);

    // 2. Pull to stencil
    stencilFld.setSize(cellCells_.size());

    forAll(cellCells_, celli)
    {
        const labelList& compactCells = cellCells_[celli];

        stencilFld[celli].setSize(compactCells.size());

        forAll(compactCells, i)
        {
            stencilFld[celli][i] = flatFld[compactCells[i]];
        }
    }
}

// ************************************************************************* //

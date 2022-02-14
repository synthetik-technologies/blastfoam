/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "locationMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::locationMapper::interpolate(Field<Type>& pf) const
{
    if (!constructMap_)
    {
        return;
    }

    Field<Type> pfOld(pf);
    pf.resize(mesh_.points().size());

    if (edgesOldPtr_.valid())
    {
        const edgeList& edges = edgesOldPtr_();
        const scalarListList& weights = edgeWeights();
        forAll(weights, i)
        {
            const edge& e = edges[edgeSplits_[i]];
            point& pt = pf[newEdgeIndices_[i]];
            pt = Zero;
            forAll(weights[i], j)
            {
                pt += pfOld[e[j]]*weights[i][j];
            }
        }
    }

    if (facesOldPtr_.valid())
    {
        const faceList& faces = facesOldPtr_();
        const scalarListList& weights = faceWeights();
        forAll(weights, i)
        {
            const face& f = faces[faceSplits_[i]];
            point& pt = pf[newFaceIndices_[i]];
            pt = Zero;
            forAll(weights[i], j)
            {
                pt += pfOld[f[j]]*weights[i][j];
            }
        }

        if (cellsOldPtr_.valid())
        {
            const cellList& cells = cellsOldPtr_();
            const scalarListList& weights = cellWeights();
            forAll(weights, i)
            {
                const cell& c = cells[cellSplits_[i]];
                const labelList cp(c.labels(faces));
                point& pt = pf[newCellIndices_[i]];
                pt = Zero;
                forAll(weights[i], j)
                {
                    pt += pfOld[cp[j]]*weights[i][j];
                }
            }
        }
    }
}

// ************************************************************************* //

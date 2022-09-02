/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "voxelMeshSearch.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Container, class Type>
void Foam::voxelMeshSearch::fill
(
    Container& elems,
    const boundBox& bb,
    const labelVector& nDivs,
    const boundBox& subBb,
    const Type val
)
{
    labelVector minIds(index3(bb, nDivs, subBb.min()));
    labelVector maxIds(index3(bb, nDivs, subBb.max()));

    for (direction cmpt = 0; cmpt < 3; cmpt++)
    {
        if (maxIds[cmpt] < 0 || minIds[cmpt] >= nDivs[cmpt])
        {
            return;
        }
        // Clip
        maxIds[cmpt] = min(maxIds[cmpt], nDivs[cmpt]-1);
        minIds[cmpt] = max(minIds[cmpt], 0);
    }

    const labelVector off(offset(nDivs));
    label voxeli = index(nDivs, minIds);
    for (label k = minIds[2]; k <= maxIds[2]; k++)
    {
        const label start1 = voxeli;
        for (label j = minIds[1]; j <= maxIds[1]; j++)
        {
            const label start0 = voxeli;
            for (label i = minIds[0]; i <= maxIds[0]; i++)
            {
                elems[voxeli] = val;
                voxeli += off[0];
            }
            voxeli = start0 + off[1];
        }
        voxeli = start1 + off[2];
    }
}


template<class Container, class Type, class CombineOp>
void Foam::voxelMeshSearch::fill
(
    Container& elems,
    const boundBox& bb,
    const labelVector& nDivs,
    const boundBox& subBb,
    const Type val,
    const CombineOp& cop
)
{
    labelVector minIds(index3(bb, nDivs, subBb.min()));
    labelVector maxIds(index3(bb, nDivs, subBb.max()));

    for (direction cmpt = 0; cmpt < 3; cmpt++)
    {
        if (maxIds[cmpt] < 0 || minIds[cmpt] >= nDivs[cmpt])
        {
            return;
        }
        // Clip
        maxIds[cmpt] = min(maxIds[cmpt], nDivs[cmpt]-1);
        minIds[cmpt] = max(minIds[cmpt], 0);
    }

    const labelVector off(offset(nDivs));
    label voxeli = index(nDivs, minIds);
    for (label k = minIds[2]; k <= maxIds[2]; k++)
    {
        const label start1 = voxeli;
        for (label j = minIds[1]; j <= maxIds[1]; j++)
        {
            const label start0 = voxeli;
            for (label i = minIds[0]; i <= maxIds[0]; i++)
            {
                cop(elems[voxeli], val);
                voxeli += off[0];
            }
            voxeli = start0 + off[1];
        }
        voxeli = start1 + off[2];
    }
}


template<class Container, class Type>
bool Foam::voxelMeshSearch::overlaps
(
    const boundBox& bb,
    const labelVector& nDivs,
    const boundBox& subBb,
    const Container& elems,
    const Type val,
    const bool isNot
)
{
    // Checks if subBb overlaps any voxel set to val

    labelVector minIds(index3(bb, nDivs, subBb.min()));
    labelVector maxIds(index3(bb, nDivs, subBb.max()));

    for (direction cmpt = 0; cmpt < 3; cmpt++)
    {
        if (maxIds[cmpt] < 0 || minIds[cmpt] >= nDivs[cmpt])
        {
            return false;
        }
        // Clip
        maxIds[cmpt] = min(maxIds[cmpt], nDivs[cmpt]-1);
        minIds[cmpt] = max(minIds[cmpt], 0);
    }

    if (elems.size() != cmptProduct(nDivs))
    {
        FatalErrorInFunction
            << "sizes:" << elems.size() << " and " << nDivs
            << exit(FatalError);
    }


    const labelVector off(offset(nDivs));
    label voxeli = index(nDivs, minIds);
    for (label k = minIds[2]; k <= maxIds[2]; k++)
    {
        const label start1 = voxeli;
        for (label j = minIds[1]; j <= maxIds[1]; j++)
        {
            const label start0 = voxeli;
            for (label i = minIds[0]; i <= maxIds[0]; i++)
            {
                const Type elemVal = elems[voxeli];
                if (isNot != (elemVal == val))
                {
                    return true;
                }
                voxeli += off[0];
            }
            voxeli = start0 + off[1];
        }
        voxeli = start1 + off[2];
    }
    return false;
}


template<class Container, class Type>
void Foam::voxelMeshSearch::write
(
    OBJstream& os,
    const boundBox& bb,
    const labelVector& nDivs,
    const Container& elems,
    const Type val,
    const bool isNot
)
{
    if (elems.size() != cmptProduct(nDivs))
    {
        FatalErrorInFunction
            << "sizes:" << elems.size() << " and " << nDivs
            << exit(FatalError);
    }

    const labelVector off(offset(nDivs));
    label voxeli = index(nDivs, labelVector(0, 0, 0));
    for (label k = 0; k < nDivs[2]; k++)
    {
        const label start1 = voxeli;
        for (label j = 0; j < nDivs[1]; j++)
        {
            const label start0 = voxeli;
            for (label i = 0; i < nDivs[0]; i++)
            {
                const Type& elemVal = elems[voxeli];
                if (isNot != (elemVal == val))
                {
                    os.write(centre(bb, nDivs, labelVector(i, j, k)));
                }
                voxeli += off[0];
            }
            voxeli = start0 + off[1];
        }
        voxeli = start1 + off[2];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

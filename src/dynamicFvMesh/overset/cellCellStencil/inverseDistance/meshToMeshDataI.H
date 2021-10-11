/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::meshToMeshData::meshToMeshData()
:
    tgtCelli_(-2)
{}


inline Foam::meshToMeshData::meshToMeshData(const label tgtCelli)
:
    tgtCelli_(tgtCelli)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::meshToMeshData::valid(TrackingData& td) const
{
    return tgtCelli_ != -2;
}


// No geometric data so never any problem on cyclics
template<class TrackingData>
inline bool Foam::meshToMeshData::sameGeometry
(
    const polyMesh&,
    const meshToMeshData&,
    const scalar,
    TrackingData& td
) const
{
    return true;
}


template<class TrackingData>
inline void Foam::meshToMeshData::leaveDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFacei,
    const point& faceCentre,
    TrackingData& td
)
{}


template<class TrackingData>
inline void Foam::meshToMeshData::transform
(
    const polyMesh&,
    const tensor& rotTensor,
    TrackingData& td
)
{}


// Update absolute geometric quantities.
template<class TrackingData>
inline void Foam::meshToMeshData::enterDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFacei,
    const point& faceCentre,
    TrackingData& td
)
{}


// Update cell with neighbouring face information
template<class TrackingData>
inline bool Foam::meshToMeshData::updateCell
(
    const polyMesh& mesh,
    const label thisCelli,
    const label neighbourFacei,
    const meshToMeshData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        const point& cc = mesh.cellCentres()[thisCelli];

        // See if we are in the neighbour
        const label tgti = neighbourInfo.tgtCelli_;

        // Try match of neighbour
        if (td.tgtMesh_.pointInCell(cc, tgti, polyMesh::CELL_TETS))
        {
            tgtCelli_ = tgti;
            return true;
        }

        // Try match of neighbour's neighbours
        const labelList& tgtNbrs = td.tgtMesh_.cellCells(tgti);
        for (const label tgti : tgtNbrs)
        {
            // Try match of neighbour
            if (td.tgtMesh_.pointInCell(cc, tgti, polyMesh::CELL_TETS))
            {
                tgtCelli_ = tgti;
                return true;
            }
        }

        // Leave unmatched cell for visit from different face
    }

    return false;
}


// Update face with neighbouring cell information
template<class TrackingData>
inline bool Foam::meshToMeshData::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label neighbourCelli,
    const meshToMeshData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        tgtCelli_ = neighbourInfo.tgtCelli_;
        return true;
    }

    return false;
}


// Update face with coupled face information
template<class TrackingData>
inline bool Foam::meshToMeshData::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const meshToMeshData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        tgtCelli_ = neighbourInfo.tgtCelli_;
        return true;
    }

    return false;
}


template<class TrackingData>
inline bool Foam::meshToMeshData::equal
(
    const meshToMeshData& rhs,
    TrackingData& td
) const
{
    if (!valid(td))
    {
        return (!rhs.valid(td));
    }
    else
    {
        return operator==(rhs);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::meshToMeshData::operator==
(
    const meshToMeshData& rhs
) const
{
    return tgtCelli_ == rhs.tgtCelli_;
}


inline bool Foam::meshToMeshData::operator!=
(
    const meshToMeshData& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //

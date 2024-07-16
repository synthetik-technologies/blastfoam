/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "feBoundaryMesh.H"
#include "lduSchedule.H"
#include "feMesh.H"
#include "globalMeshData.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::feBoundaryMesh::feBoundaryMesh
(
    const feMesh& m,
    const polyBoundaryMesh& basicBdry
)
:
    fePatchList(basicBdry.size()),
    mesh_(m)
{
    reset(basicBdry);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::feBoundaryMesh::reset(const polyBoundaryMesh& basicBdry)
{
    // Set boundary patches
    fePatchList& patches = *this;

    forAll(patches, patchi)
    {
        patches.set
        (
            patchi,
            fePatch::New(basicBdry[patchi], *this).ptr()
        );
    }
}



void Foam::feBoundaryMesh::calcGeometry()
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initCalcGeometry(pBufs);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).calcGeometry(pBufs);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        forAll(patchSchedule, patchEvali)
        {
            label patchi = patchSchedule[patchEvali].patch;

            if (patchSchedule[patchEvali].init)
            {
                operator[](patchi).initCalcGeometry(pBufs);
            }
            else
            {
                operator[](patchi).calcGeometry(pBufs);
            }
        }
    }
}


Foam::label Foam::feBoundaryMesh::findPatchID(const word& patchName) const
{
   return mesh().mesh().boundaryMesh().findPatchID(patchName);
}


Foam::labelList Foam::feBoundaryMesh::findIndices
(
    const wordRe& key,
    const bool usePatchGroups
) const
{
    return mesh().mesh().boundaryMesh().findIndices(key, usePatchGroups);
}


void Foam::feBoundaryMesh::movePoints(const pointField& p)
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initMovePoints(pBufs, p);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).movePoints(pBufs, p);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        forAll(patchSchedule, patchEvali)
        {
            label patchi = patchSchedule[patchEvali].patch;

            if (patchSchedule[patchEvali].init)
            {
                operator[](patchi).initMovePoints(pBufs, p);
            }
            else
            {
                operator[](patchi).movePoints(pBufs, p);
            }
        }
    }
}


void Foam::feBoundaryMesh::shuffle
(
    const labelUList& newToOld,
    const bool validBoundary
)
{
    fePatchList& patches = *this;
    patches.shuffle(newToOld);
}



// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

const Foam::fePatch& Foam::feBoundaryMesh::operator[]
(
    const word& patchName
) const
{
    const label patchi = findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch named " << patchName << " not found." << nl
            << abort(FatalError);
    }

    return operator[](patchi);
}


Foam::fePatch& Foam::feBoundaryMesh::operator[]
(
    const word& patchName
)
{
    const label patchi = findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch named " << patchName << " not found." << nl
            << abort(FatalError);
    }

    return operator[](patchi);
}


// ************************************************************************* //

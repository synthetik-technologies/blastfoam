/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenCFD Ltd.
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

Class
    Foam::fvMeshPrimitiveLduAddressing

Description
    Variant of fvMeshLduAddressing that contains addressing instead of
    slices.

    Alternatively use lduPrimitiveMesh but that assumes there are only
    lduInterfaces and not also generic patches.

SourceFiles
    fvMeshPrimitiveLduAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshPrimitiveLduAddressing_H
#define fvMeshPrimitiveLduAddressing_H

#include "lduPrimitiveMesh.H"
#include "lduAddressing.H"
#include "fvMesh.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class fvMeshPrimitiveLduAddressing Declaration
\*---------------------------------------------------------------------------*/

class fvMeshPrimitiveLduAddressing
:
    public lduAddressing
{
    // Private data

        //- Lower (face to owner addressing)
        labelList lowerAddr_;

        //- Upper (face to neighbour addressing)
        labelList upperAddr_;

        //- Patch addressing as a list of sublists
        List<const labelUList*> patchAddr_;

        //- Patch field evaluation schedule
        const lduSchedule& patchSchedule_;


    // Private Member Functions

        static labelList upperTriOrder
        (
            const label nCells,
            const labelUList& lower,
            const labelUList& upper
        );

        //- No copy construct
        fvMeshPrimitiveLduAddressing
        (
            const fvMeshPrimitiveLduAddressing&
        ) = delete;

        //- No copy assignment
        void operator=(const fvMeshPrimitiveLduAddressing&) = delete;


public:

    // Constructors

        //- Construct from mesh
        fvMeshPrimitiveLduAddressing(const fvMesh& mesh);

        //- Construct from components
        fvMeshPrimitiveLduAddressing
        (
            const label nCells,
            labelList&& lowerAddr,
            labelList&& upperAddr,
            const List<const labelUList*>& interfaces,
            const lduSchedule& ps
        );


    //- Destructor
    virtual ~fvMeshPrimitiveLduAddressing() = default;


    // Member Functions

        //- Return lower addressing (i.e. lower label = upper triangle)
        virtual const labelUList& lowerAddr() const
        {
            return lowerAddr_;
        }

        //- Return upper addressing (i.e. upper label)
        virtual const labelUList& upperAddr() const
        {
            return upperAddr_;
        }

        //- Return patch addressing
        virtual const labelUList& patchAddr(const label i) const
        {
            return *patchAddr_[i];
        }

        // Return patch field evaluation schedule
        virtual const lduSchedule& patchSchedule() const
        {
            return patchSchedule_;
        }

        //- Given additional addressing (in the form of additional neighbour
        //  cells, i.e. like cellCells)
        //  - add any additional faces
        //  - sort in upper-triangular order
        //  - construct cell-faces equivalent of given nbrCells
        //    (so e.g. nbrCellFaces[cellI][0] is the face corresponding to
        //     the cell nbrCells[cellI][0])
        //    (note: face in nbrCellFaces is -1 if the nbrCell is not local)
        //  - construct additional processor interface addressing:
        //      per processor the local and the remote cell.
        //  - return old-to-new face mapping
        static labelList addAddressing
        (
            const lduAddressing& addr,
            const labelListList& nbrCells,
            label& nExtraFaces,
            labelList& lower,
            labelList& upper,
            labelListList& nbrCellFaces,
            const globalIndex&,
            const labelList& globalCellIDs,
            labelListList& localFaceCells,
            labelListList& remoteFaceCells
        );

        //- Return off-diagonal index given owner and neighbour label. Return
        //  -1 if not found
        static label triIndex(const lduAddressing&, const label, const label);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

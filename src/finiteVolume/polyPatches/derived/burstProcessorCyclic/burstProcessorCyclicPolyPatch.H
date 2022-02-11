/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Class
    Foam::burstProcessorCyclicPolyPatch

Description
    Poly patch for a burstProcessorCyclicing surface (cannot be undone)

SourceFiles
    burstProcessorCyclicPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef burstProcessorCyclicPolyPatch_H
#define burstProcessorCyclicPolyPatch_H

#include "processorCyclicPolyPatch.H"
#include "burstCyclicPolyPatch.H"
#include "burstPolyPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class burstProcessorCyclicPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class burstProcessorCyclicPolyPatch
:
    public processorCyclicPolyPatch,
    public burstPolyPatchBase
{
public:

    //- Runtime type information
    TypeName("burstProcessorCyclic");


    // Constructors

        //- Construct from components
        burstProcessorCyclicPolyPatch
        (
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const int myProcNo,
            const int neighbProcNo,
            const word& referPatchName,
            const word& patchType = typeName
        );

        //- Construct from dictionary
        burstProcessorCyclicPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        burstProcessorCyclicPolyPatch
        (
            const burstProcessorCyclicPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        burstProcessorCyclicPolyPatch
        (
            const burstProcessorCyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct as given the original patch and resetting the
        //  face list, boundary mesh information and referPatch
        burstProcessorCyclicPolyPatch
        (
            const burstProcessorCyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& referPatchName
        );

        //- Construct given the original patch and a map
        burstProcessorCyclicPolyPatch
        (
            const burstProcessorCyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new burstProcessorCyclicPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new burstProcessorCyclicPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& referPatchName
        ) const
        {
            return autoPtr<polyPatch>
            (
                new burstProcessorCyclicPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    referPatchName
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new burstProcessorCyclicPolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~burstProcessorCyclicPolyPatch();


    // Member Functions

        // Access

            const burstCyclicPolyPatch& referPatch() const
            {
                const polyPatch& pp = this->boundaryMesh()[referPatchID()];
                return refCast<const burstCyclicPolyPatch>(pp);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

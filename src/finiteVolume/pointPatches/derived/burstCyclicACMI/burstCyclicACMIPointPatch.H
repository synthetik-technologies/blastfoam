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
    Foam::burstCyclicACMIPointPatch

Description
    Cyclic patch for post-processing.

SourceFiles
    burstCyclicACMIPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef burstCyclicACMIPointPatch_H
#define burstCyclicACMIPointPatch_H

#include "cyclicACMIPointPatch.H"
#include "burstCyclicACMIPolyPatch.H"
#include "pointBoundaryMesh.H"
#include "valuePointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class burstCyclicACMIPointPatch Declaration
\*---------------------------------------------------------------------------*/

class burstCyclicACMIPointPatch
:
    public cyclicACMIPointPatch
{
    // Private Data

        //- Const reference to the burstCyclicACMIPolyPatch
        const burstCyclicACMIPolyPatch& burstCyclicACMIPolyPatch_;


    // Private Member Functions

        //- Disallow default bitwise copy construction
        burstCyclicACMIPointPatch(const burstCyclicACMIPointPatch&);

        //- Disallow default bitwise assignment
        void operator=(const burstCyclicACMIPointPatch&);


    // Demand driven private data

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName(burstCyclicACMIPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        burstCyclicACMIPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    //- Destructor
    virtual ~burstCyclicACMIPointPatch();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

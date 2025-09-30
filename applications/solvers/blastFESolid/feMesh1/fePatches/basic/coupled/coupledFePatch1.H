/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
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

Class
    Foam::coupledFePatch1

Description
    coupled patch for post-processing.  Used as the base class for processor
    and cyclic pointPatches

SourceFiles
    coupledFePatch1.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFePatch1_H
#define coupledFePatch1_H

#include "fePatch1.H"
#include "coupledPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class coupledFePatch Declaration
\*---------------------------------------------------------------------------*/

class coupledFePatch1
:
    public fePatch1
{
    // Private Data

        const coupledPolyPatch& coupledPolyPatch_;


    // Private Member Functions

        //- Disallow default bitwise copy construction
        coupledFePatch1(const coupledFePatch1&);

        //- Disallow default bitwise assignment
        void operator=(const coupledFePatch1&);


protected:

    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&) = 0;

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&) = 0;

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&) = 0;

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&) = 0;

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&) = 0;

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&) = 0;


public:

    //- Runtime type information
    TypeName(coupledPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        coupledFePatch1
        (
            const polyPatch& patch,
            const feBoundaryMesh1& bm
        );


    //- Destructor
    virtual ~coupledFePatch1();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

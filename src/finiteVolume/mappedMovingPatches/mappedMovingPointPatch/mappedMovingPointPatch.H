/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-08-21 Synthetik Applied Technologies: Mapping of point patches
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
    Foam::mappedMovingPointPatch

Description
    Determines a mapping between patch points and processors they're
    on.

SourceFiles
    mappedMovingPointPatch.C
    mappedMovingPointPatchBase.H
    mappedMovingPointPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef mappedMovingPointPatch_H
#define mappedMovingPointPatch_H

#include "facePointPatch.H"
#include "mappedMovingPointPatchBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointMesh;

/*---------------------------------------------------------------------------*\
                      Class mappedMovingPointPatch Declaration
\*---------------------------------------------------------------------------*/

class mappedMovingPointPatch
:
    public facePointPatch,
    public mappedMovingPointPatchBase
{

protected:

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
    TypeName("mappedMovingWall");


    // Constructors

        //- Construct from components
        mappedMovingPointPatch
        (
            const polyPatch& pp,
            const pointBoundaryMesh& bm
        );

    //- Destructor
    virtual ~mappedMovingPointPatch();


    // Member Functions

        //- Write the PointPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

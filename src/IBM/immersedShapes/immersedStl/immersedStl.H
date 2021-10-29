/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technology
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::immersedStl

Description


SourceFiles
    immersedStl.C

\*---------------------------------------------------------------------------*/

#ifndef immersedStl_H
#define immersedStl_H

#include "immersedShape.H"
#include "Time.H"
#include "triSurfaceMesh.H"
#include "edgeIntersections.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                       Class immersedStl Declaration
\*---------------------------------------------------------------------------*/

class immersedStl
:
    public immersedShape
{
    // Private data

        //- Save the constructing dictionary
        const dictionary dict_;

        //- Name of stl file
        word fileName_;

        //- Immersed boundary trisurface mesh
        triSurfaceMesh immersedMesh_;

        //- Class used to search original triMesh
        autoPtr<triSurfaceSearch> tssPtr_;


    // Private member functions

        //- Coarsen mesh given a upper and lower values of dx (2-D)
        void coarsenRefine2D(pointField& points);

        //- Return 2D points
        void get2DPoints(List<face>& faces);

        // Functions for calculating intersections with mesh
        bool intersectSurfaces
        (
            triSurface& surf,
            edgeIntersections& edgeCuts
        );

        bool intersectSurfaces
        (
            triSurface& surf1,
            edgeIntersections& edgeCuts1,
            triSurface& surf2,
            edgeIntersections& edgeCuts2
        );

        void calcEdgeCuts
        (
            triSurface& surf1,
            triSurface& surf2,
            const bool perturb,
            edgeIntersections& edge1Cuts,
            edgeIntersections& edge2Cuts
        );

        void addAxisPoints(pointField& points) const;

        void createRegionTriSurfaces
        (
            const triSurface& triMesh,
            PtrList<triSurface>& regionTriSurfaces,
            boolList& onAxis
        ) const;

        //- Does the centre of mass need to be read
        virtual bool needCentreOfMass() const
        {
            return true;
        }

        //- Does the orientation need to be read
        virtual bool needOrientation() const
        {
            return true;
        }

public:

    //- Runtime type information
    TypeName("stl");


    // Constructors

        //- Construct from patch
        immersedStl
        (
            const polyMesh&,
            const immersedBoundaryObject&,
            const dictionary& dict
        );



    //- Destructor
    virtual ~immersedStl();


    // Member Functions

        //- Return bool list if list of points is inside the object
        virtual labelList calcInside(const pointField& points) const;

        //- Return if point is inside the object
        virtual bool inside(const point& pt) const;

        //- Return moment of inertia
        virtual diagTensor momentOfInertia() const
        {
            if (ei_ == -1)
            {
                return dict_.lookup<diagTensor>("momentOfInertia");
            }
            diagTensor dt(great, great, great);
            dt[ei_] = diagTensor
            (
                dict_.lookup<diagTensor>("momentOfInertia")
            )[ei_];
            return dt;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::immersedFvMesh

Description


SourceFiles
    immersedFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef immersedFvMesh_H
#define immersedFvMesh_H

#include "immersedShape.H"
#include "triSurfaceMesh.H"
#include "dynamicBlastFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class immersedFvMesh Declaration
\*---------------------------------------------------------------------------*/

class immersedFvMesh
:
    public immersedShape
{
    // Private data

        //- Immersed fvMesh
        autoPtr<dynamicBlastFvMesh> fvMeshPtr_;
        dynamicBlastFvMesh& fvMesh_;

        //- Patch name
        word patchName_;

        //- Full (non parallel) triangulated mesh that includes all patches
        //  except for processor patches
        //  Used to check if a point is internal
        triSurface* fullTriMeshPtr_;

        autoPtr<triSurfaceSearch> tssPtr_;


    // Private member functions

        //- Coarsen mesh given a upper and lower values of dx (2-D)
//         void coarsenRefine2D(pointField& points);

        //- Return 2D points
//         void get2DPoints(List<labelledTri>& faces);

        //- Does the centre of mass need to be read
        virtual bool needCentreOfMass() const
        {
            return false;
        }

        //- Triangulate the mesh given patches
        autoPtr<triSurface> triangulate
        (
            const polyBoundaryMesh& bMesh,
            const labelHashSet& patches
        ) const;

        //- Does the orientation need to be read
        virtual bool needOrientation() const
        {
            return false;
        }

public:

    //- Runtime type information
    TypeName("immersedFvMesh");


    // Constructors

        //- Construct from patch
        immersedFvMesh
        (
            const polyMesh&,
            const immersedBoundaryObject&,
            const dictionary& dict
        );



    //- Destructor
    virtual ~immersedFvMesh();


    // Member Functions

        //- Return the fvMesh
        const fvMesh& immersedMesh() const
        {
            return fvMesh_;
        }

        //- Access the fvMesh
        fvMesh& immersedMesh()
        {
            return fvMesh_;
        }

        //- Return the fvMesh
        const dynamicBlastFvMesh& immersedDyMesh() const
        {
            return fvMesh_;
        }

        //- Access the fvMesh
        dynamicBlastFvMesh& immersedDyMesh()
        {
            return fvMesh_;
        }

        //- Move the triMesh
        virtual void movePoints();

        //- Return bool list if list of points is inside the object
        virtual labelList calcInside(const pointField& points) const;

        //- Return if point is inside the object
        virtual bool inside(const point& pt) const;

        //- Return moment of inertia
        virtual diagTensor momentOfInertia() const
        {
//             if (ei_ == -1)
//             {
//                 return dict_.lookup<diagTensor>("momentOfInertia");
//             }
            diagTensor dt(great, great, great);
//             dt[ei_] = diagTensor
//             (
//                 dict_.lookup<diagTensor>("momentOfInertia")
//             )[ei_];
            return dt;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

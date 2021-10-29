/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
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
    Foam::immersedShape

Description


SourceFiles
    immersedShape.C

\*---------------------------------------------------------------------------*/

#ifndef immersedShape_H
#define immersedShape_H

#include "pointField.H"
#include "polyMesh.H"
#include "triSurface.H"
#include "standAlonePatch.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class immersedBoundaryObject;

/*---------------------------------------------------------------------------*\
                       Class immersedShape Declaration
\*---------------------------------------------------------------------------*/

class immersedShape
{
protected:
// protected:

    // Protected data

        //- Reference to mesh
        const polyMesh& pMesh_;

        //- Reference to the immersed boundary object
        const immersedBoundaryObject& object_;

        //- Desired mesh size
        scalar dx_;

        //- Points
        pointField points0_;

        //- Old face centres
        pointField faceCentresOld_;

        //- Used geometric dimensions
        vector geometricD_;

        // Axisymmetric Variables

            //- Radial index
            label ri_;

            //- Axis index
            label ai_;

            //- Angle of wedge
            scalar angle_;

            //- Scale factor used to scale mass
            scalar scale_;


        //- Empty geometric directions index
        label ei_;

        //- "x" directions
        label xi_;

        //- "y" direction
        label yi_;

        //- Centre of mass
        vector centreOfMass_;

        //- Orientation of shape
        tensor orientation_;

        //- Name of shape
        word name_;

        //- Write sample when mesh has been created
        Switch write_;

        //- Immersed boundary trisurface mesh created from shape
        autoPtr<standAlonePatch> patchPtr_;

        //- Bounding boundBox
        boundBox bb_;



    // Protected Member Functions

        //- Return rotation tensor
        tensor RotM(const scalar) const;


        // 2-D and Axisymmetric functions

            //- Sort points counter-clockwise
            tmp<scalarField> sortPointsPolar(pointField& points) const;

            //- Sort points counter-clockwise
            void sortPointsConcave(pointField& points);

            //- Return points inside bound box
            void insidePoints
            (
                pointField& points,
                labelList& indices
            ) const;

            //- Make list of points axisymmetric and create face indexing
            void extrudeAxi
            (
                pointField& points,
                faceList& faces,
                const bool onAxis = true
            ) const;

            //- Extrude Points (2-D)
            void extrude2D
            (
                pointField& points,
                faceList& faces
            ) const;

            //- Refine a 3-D triMesh
            void refine3D(triSurface& triMesh);

            //- Adjust points
            virtual bool adjustPoints(pointField& points)
            {
                return false;
            }

            //- Does the centre of mass need to be read
            virtual bool needCentreOfMass() const = 0;

            //- Does the orientation need to be read
            virtual bool needOrientation() const = 0;

        //- Read centre of mass and orientation if needed
        void read(const dictionary&);


public:

    //- Runtime type information
    TypeName("immersedShape");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            immersedShape,
            twoD,
            (
                const polyMesh& mesh,
                const immersedBoundaryObject& ibo,
                const dictionary& dict
            ),
            (mesh, ibo, dict)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            immersedShape,
            threeD,
            (
                const polyMesh& mesh,
                const immersedBoundaryObject& ibo,
                const dictionary& dict
            ),
            (mesh, ibo, dict)
        );


    /// Selectors

        static autoPtr<immersedShape> New
        (
            const polyMesh& mesh,
            const immersedBoundaryObject& ibo,
            const dictionary& dict
        );

        immersedShape
        (
            const polyMesh& mesh,
            const immersedBoundaryObject& ibo,
            const dictionary& dict
        );


    //- Check shape type and return non-constant access to that type
    template<class ShapeType>
    static ShapeType& castShapeType(const immersedShape& shape)
    {
        if (!isA<ShapeType>(shape))
        {
            FatalErrorInFunction
                << "Invalid immersed shape specified."
                << "Expected " << ShapeType::typeName
                << " but got " << shape.typeName
                << endl
                << abort(FatalError);
        }
        return const_cast<ShapeType&>(refCast<const ShapeType>(shape));
    }


    //- Destructor
    virtual ~immersedShape();


    // Member Functions

        // Access

            //- Return name
            const word& name() const
            {
                return name_;
            }

            //- Return the scale factor
            const scalar& scale() const
            {
                return scale_;
            }

            //- Return points
            const pointField& points() const
            {
                return patchPtr_->localPoints();
            }

            //- Return old points
            const pointField& faceCentresOld() const
            {
                return faceCentresOld_;
            }

            //- Return moment of inertia
            virtual diagTensor momentOfInertia() const = 0;

            //- Return the triSurface
            const standAlonePatch& patch() const
            {
                return patchPtr_();
            }

            //- Access to the triSurface
            standAlonePatch& patch()
            {
                return patchPtr_();
            }

            //- Return the bounding box
            const boundBox& bounds() const
            {
                return bb_;
            }

            //- Return the shapes centre of mass
            virtual vector centreOfMass() const
            {
                return centreOfMass_;
            }

            //- Return the initial orientation of the shape
            virtual tensor orientation() const
            {
                return orientation_;
            }


        //- Return bool list if list of points is inside the object
        virtual labelList calcInside(const pointField& points) const = 0;

        //- Return if point is inside the object
        virtual bool inside(const point& pt) const = 0;

        //- Move the triMesh
        virtual void movePoints();

        //- Zero empty direction
        vector zeroDir(const vector& v) const
        {
            vector vNew(v);
            if (ei_ != -1)
            {
                vNew[ei_] = 0;
            }
            return vNew;
        }

        //- Zero empty direction
        void addMagVector(const vector& v, scalar& magVTot, label& n) const
        {
            if (ei_ != -1)
            {
                scalar magX(mag(v[xi_]));
                scalar magY(mag(v[yi_]));
                scalar magE(mag(v[ei_]));
                if (magE < magX || magE < magY)
                {
                    magVTot += mag(v);
                    n++;
                }
                return;
            }
            magVTot += mag(v);
            n++;
        }

        //- Calculate line slope and y intercept given two points
        //  Return true if the line is vertical
        bool calcLine
        (
            const vector& p1,
            const vector& p2,
            scalar& m,
            scalar& b
        ) const;

        //- Do lines intersect ([p1, p2], [p3, p4])
        bool intersection
        (
            const vector& p1, const vector& p2,
            const vector& p3, const vector& p4
        ) const;

        //- Calculate line and normal
        vector calcNormal
        (
            const vector& p1,
            const vector& p2
        ) const;

        //- Make sure normals are pointing towards a given point
        void correctNormal
        (
            vector& normal,
            const point& start,
            const point& end,
            const point& p1,
            const point& p2
        ) const;

        //- Create list of points given a start and end point
        tmp<pointField> discretizeLine
        (
            const vector& p1,
            const vector& p2,
            const bool addEnd = false
        ) const;

        //-  Write the mesh VTK
        void writeVTK() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::immersedEllipsoid

Description


SourceFiles
    immersedEllipsoid.C

\*---------------------------------------------------------------------------*/

#ifndef immersedEllipsoid_H
#define immersedEllipsoid_H

#include "immersedShape.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                       Class immersedEllipsoid Declaration
\*---------------------------------------------------------------------------*/

class immersedEllipsoid
:
    public immersedShape
{
    // Private data

        //- First axis radius
        scalar a_;

        //- Second axis radius
        scalar b_;

        //- Third axis radius
        scalar c_;

        //- Move points to the surface of the sphere
        virtual bool adjustPoints(pointField& points);

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
    TypeName("ellipsoid");


    // Constructors

        //- Construct from patch
        immersedEllipsoid
        (
            const polyMesh&,
            const immersedBoundaryObject&,
            const dictionary& dict
        );



    //- Destructor
    virtual ~immersedEllipsoid();


    // Member Functions

        //- Return bool list if list of points is inside the object
        virtual labelList calcInside(const pointField& points) const;

        //- Return if point is inside the object
        virtual bool inside(const point& pt) const;

        //- Return moment of inertia
        virtual diagTensor momentOfInertia() const
        {
            return diagTensor
            (
                1.0/5.0*(sqr(b_) + sqr(c_)),
                1.0/5.0*(sqr(a_) + sqr(c_)),
                1.0/5.0*(sqr(a_) + sqr(b_))
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

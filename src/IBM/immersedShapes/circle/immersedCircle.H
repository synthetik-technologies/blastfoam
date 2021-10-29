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
    Foam::immersedCircle

Description


SourceFiles
    immersedCircle.C

\*---------------------------------------------------------------------------*/

#ifndef immersedCircle_H
#define immersedCircle_H

#include "immersedShape.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                       Class immersedCircle Declaration
\*---------------------------------------------------------------------------*/

class immersedCircle
:
    public immersedShape
{
    // Private data

        //- Radius
        scalar radius_;

        //- Return 2D points
        virtual void get2DPoints();

        //- Does the centre of mass need to be read
        virtual bool needCentreOfMass() const
        {
            return true;
        }

        //- Does the orientation need to be read
        virtual bool needOrientation() const
        {
            return false;
        }


public:

    //- Runtime type information
    TypeName("circle");


    // Constructors

        //- Construct from patch
        immersedCircle
        (
            const polyMesh&,
            const immersedBoundaryObject&,
            const dictionary& dict
        );



    //- Destructor
    virtual ~immersedCircle();


    // Member Functions

        //- Return bool list if list of points is inside the object
        virtual labelList calcInside(const pointField& points) const;

        //- Return if point is inside the object
        virtual bool inside(const point& pt) const;

        //- Return moment of inertia
        virtual diagTensor momentOfInertia() const
        {
            scalar pi = Foam::constant::mathematical::pi;
            diagTensor dt(great, great, great);
            dt[xi_] = 0.25*pi*pow4(radius_);
            dt[yi_] = 0.25*pi*pow4(radius_);
            dt[ei_] = 0.5*pi*pow4(radius_);
            return dt;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

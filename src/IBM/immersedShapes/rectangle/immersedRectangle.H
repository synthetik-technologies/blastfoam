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
    Foam::immersedRectangle

Description


SourceFiles
    immersedRectangle.C

\*---------------------------------------------------------------------------*/

#ifndef immersedRectangle_H
#define immersedRectangle_H

#include "immersedShape.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                       Class immersedRectangle Declaration
\*---------------------------------------------------------------------------*/

class immersedRectangle
:
    public immersedShape
{
    // Private data

        //- Length
        scalar L_;

        //- Height
        scalar H_;

        //- List of points
        pointField points_;

        //- Bounding box of unmoved rectangle
        boundBox origBb_;


public:

    //- Runtime type information
    TypeName("rectangle");


    // Constructors

        //- Construct from patch
        immersedRectangle
        (
            const polyMesh&,
            const immersedBoundaryObject&,
            const dictionary& dict
        );



    //- Destructor
    virtual ~immersedRectangle();


    // Member Functions

        //- Return a patch created with the given geometry
        virtual autoPtr<standAlonePatch> createPatch() const;

        //- Return bool list if list of points is inside the object
        virtual labelList calcInside(const pointField& points) const;

        //- Return if point is inside the object
        virtual bool inside(const point& pt) const;

        //- Return moment of inertia
        virtual diagTensor momentOfInertia() const
        {
            diagTensor dt(great, great, great);
            dt[xi_] = 1.0/12.0*L_*pow3(H_);
            dt[yi_] = 1.0/12.0*H_*pow3(L_);
            dt[ei_] = 1.0/12.0*L_*H_*(sqr(L_) + sqr(H_));
            return dt;
        }

        // Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

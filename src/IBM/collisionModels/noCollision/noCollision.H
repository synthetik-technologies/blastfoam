/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::noCollision

Description
    Place holder for 'none' option

SourceFiles
    noCollision.C

\*---------------------------------------------------------------------------*/

#ifndef immersedNoCollision_H
#define immersedNoCollision_H

#include "collisionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class noCollision Declaration
\*---------------------------------------------------------------------------*/

class noCollision
:
    public collisionModel
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noCollision
        (
            const polyMesh& mesh,
            PtrList<immersedBoundaryObject>& objects,
            const dictionary& dict
        );


    //- Destructor
    virtual ~noCollision();


    // Member Functions

        //- Collision function
        virtual void collide();

        //- Collide with walls
        virtual void updateExternalForce();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

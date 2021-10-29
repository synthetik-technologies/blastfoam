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
    Foam::collisionModel

Description
    Templated collision model class.

SourceFiles
    collisionModel.C
    collisionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef immersedCollisionModel_H
#define immersedCollisionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "immersedBoundaryObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class collisionModel Declaration
\*---------------------------------------------------------------------------*/

class collisionModel
{
protected:
    // Protected data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- List of immersed object
        PtrList<immersedBoundaryObject>& objects_;


public:

    //- Runtime type information
    TypeName("collisionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        collisionModel,
        dictionary,
        (
            const polyMesh& mesh,
            PtrList<immersedBoundaryObject>& objects,
            const dictionary& dict
        ),
        (mesh, objects, dict)
    );


    // Constructors

        //- Construct from components
        collisionModel
        (
            const polyMesh& mesh,
            PtrList<immersedBoundaryObject>& objects,
            const dictionary& dict
        );


    //- Destructor
    virtual ~collisionModel();


    //- Selector
    static autoPtr<collisionModel> New
    (
        const polyMesh& mesh,
        PtrList<immersedBoundaryObject>& objects,
        const dictionary& dict
    );


    // Member Functions

        //- Collision function
        virtual void collide() = 0;

        //- Collide with walls
        virtual void updateExternalForce() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

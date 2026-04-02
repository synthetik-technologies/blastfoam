/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020
     \\/     M anipulation  | Synthetik Applied Technology
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
    Foam::ownNeiSurfaceFields

Description
    Holder for own (1.0) and nei (-1.0) surface fields

SourceFiles
    ownNeiSurfaceFields.C

\*---------------------------------------------------------------------------*/

#ifndef ownNeiSurfaceFields_H
#define ownNeiSurfaceFields_H

#include "DemandDrivenMeshObject.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ownNeiSurfaceFields Declaration
\*---------------------------------------------------------------------------*/

class ownNeiSurfaceFields
:
    public DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        ownNeiSurfaceFields
    >
{
    surfaceScalarField own_;
    surfaceScalarField nei_;


public:

    friend class DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        ownNeiSurfaceFields
    >;

    //- Runtime type information
    TypeName("ownNeiSurfaceFields");


    // Constructors

        //- Construct from mesh
        ownNeiSurfaceFields(const fvMesh&);

    //- Destructor
    virtual ~ownNeiSurfaceFields();


    // Member Functions

        //- Return the owner interpolated field
        inline const surfaceScalarField& own() const
        {
            return own_;
        }

        //- Return the neighbor interpolated field
        inline const surfaceScalarField& nei() const
        {
            return nei_;
        }

        //- Move points
        virtual bool movePoints()
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

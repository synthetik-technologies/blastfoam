/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    UnsLinearGeomSolid

Description
    Base class for unstructured incremental solid models

SourceFiles
    UnsLinearGeomSolid.C

\*---------------------------------------------------------------------------*/

#ifndef UnsLinearGeomSolid_H
#define UnsLinearGeomSolid_H

#include "LinearGeomSolid.H"
#include "unsSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class UnsLinearGeomSolid Declaration
\*---------------------------------------------------------------------------*/

template<class IncrementalModel>
class UnsLinearGeomSolid
:
    public LinearGeomSolid<IncrementalModel>
{
protected:
    // Protected data

        //- Implicit stiffness surface field
        surfaceScalarField impKf_;


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        UnsLinearGeomSolid(const UnsLinearGeomSolid&);

        //- Disallow default bitwise assignment
        void operator=(const UnsLinearGeomSolid&);


public:

    // Constructors

        //- Construct from components
        UnsLinearGeomSolid
        (
            const word& type,
            dynamicFvMesh& mesh,
            const bool isSolid = true
        );

    // Destructor

        virtual ~UnsLinearGeomSolid()
        {}


    // Member Functions

        //- Update fields
        virtual void update(const bool correctSigma = true);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UnsLinearGeomSolid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

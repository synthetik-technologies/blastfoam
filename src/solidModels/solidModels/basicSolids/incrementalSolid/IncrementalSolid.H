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
    IncrementalSolid

Description
    Templated class for incremental solid models

SourceFiles
    IncrementalSolid.C

\*---------------------------------------------------------------------------*/

#ifndef IncrementalSolid_H
#define IncrementalSolid_H

#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class linGeomSolid Declaration
\*---------------------------------------------------------------------------*/

template<class SolidModel>
class IncrementalSolid
:
    public SolidModel
{
protected:
    // Protected data

        //- Implicit stiffness
        const volScalarField impK_;

        //- Implicit stiffness surface field
        surfaceScalarField impKf_;


    // Protected member functions

        //- Update displacement fields
        virtual void updateDisplacement();


public:

    // Constructors

        //- Construct from components
        IncrementalSolid
        (
            const word& type,
            dynamicFvMesh& mesh,
            const nonLinearGeometry::nonLinearType nonLinear,
            const bool isSolid = true
        );


    // Destructor

        virtual ~IncrementalSolid()
        {}


    // Member Functions

        // Access

            //- This solid model solves for DD instead of D
            virtual bool incremental() const
            {
                return true;
            }

            //- Each solidModel must indicate whether D or DD is the
            //  primary solution variable
            virtual const volVectorField& solutionD() const
            {
                // This model solves for DD
                return this->DD();
            }

            //- Each solidModel must indicate whether gradD or gradDD
            //  is the primary solution variable
            virtual const volTensorField& solutionGradD() const
            {
                // This model solves for DD
                return this->gradDD();
            }

            //- Return the patch implicit stiffness
            virtual tmp<scalarField> impK(const fvPatch& patch) const
            {
                return impK_.boundaryField()[patch.index()];
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IncrementalSolid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

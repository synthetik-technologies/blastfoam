/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    fixedSolid

Description
    Mathematical model where linear geometry is assumed i.e. small strains and
    small rotations are assumed, where the increment of displacement is the
    primary unknown.

    The stress is calculated by the run-time selectable mechanical law.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    fixedSolid.C

\*---------------------------------------------------------------------------*/

#ifndef fixedSolid_H
#define fixedSolid_H

#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class fixedSolid Declaration
\*---------------------------------------------------------------------------*/

class fixedSolid
:
    public solidModel
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        fixedSolid(const fixedSolid&);

        //- Disallow default bitwise assignment
        void operator=(const fixedSolid&);


public:

    //- Runtime type information
    TypeName("fixed");

    // Constructors

        //- Construct from components
        fixedSolid(fvMesh& mesh);


    // Destructor

        virtual ~fixedSolid()
        {}


    // Member Functions

        //- Evolve the solid solver and solve the mathematical model
        virtual bool evolve();

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::NONE;
        }

        //- Each solidModel must indicate whether D or DD is the primary
        //  solution variable
        virtual const volVectorField& solutionD() const
        {
            NotImplemented;
            return NullObjectRef<volVectorField>();
        }

        //- Each solidModel must indicate whether gradD or gradDD
        //  is the primary solution variable
        virtual const volTensorField& solutionGradD() const
        {
            NotImplemented;
            return NullObjectRef<volTensorField>();
        }

        //- Update fields
        virtual void update(const bool correctSigma = true)
        {}

        //- Boundary normal
        virtual tmp<vectorField> nf(const fvPatch& patch) const
        {
            return patch.nf();
        }

        //- Return the patch implicit stiffness
        virtual tmp<scalarField> impK(const fvPatch& patch) const
        {
            return tmp<scalarField>(new scalarField(patch.size(), great));
        }

        // Member Functions

        // Access

            //- Return const reference to D
            virtual const volVectorField& D() const
            {
                NotImplemented;
                return NullObjectRef<volVectorField>();
            }

            //- Return const reference to U
            virtual const volVectorField& U() const
            {
                NotImplemented;
                return NullObjectRef<volVectorField>();
            }

            //- Return const reference to pointD
            virtual const pointVectorField& pointD() const
            {
                NotImplemented;
                return NullObjectRef<pointVectorField>();
            }

            //- Return const reference to sigma
            virtual const volSymmTensorField& sigma() const
            {
                NotImplemented;
                return NullObjectRef<volSymmTensorField>();
            }

            //- Return the patch stress
            virtual const symmTensorField& sigma(const fvPatch& patch) const
            {
                NotImplemented;
                return NullObjectRef<symmTensorField>();
            }


            //- Return enforceLinear switch
            virtual bool enforceLinear() const
            {
                return false;
            }

            //- Return enforceLinear switch
            virtual void enforceLinear(const bool tf)
            {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

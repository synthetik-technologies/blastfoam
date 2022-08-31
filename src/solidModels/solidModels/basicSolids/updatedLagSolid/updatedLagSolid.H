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
    updatedLagSolid

Description
    Base class for updated Lagrangian solid models

SourceFiles
    updatedLagSolid.C

\*---------------------------------------------------------------------------*/

#ifndef updatedLagSolid_H
#define updatedLagSolid_H

#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class updatedLagSolid Declaration
\*---------------------------------------------------------------------------*/

template<class IncrementalModel>
class updatedLagSolid
:
    public IncrementalModel
{
protected:
    // Protected data

        //- Total deformation gradient
        volTensorField F_;

        //- Jacobian of total deformation gradient
        volScalarField J_;

        //- Relative deformation gradient
        volTensorField relF_;

        //- Relative Jacobian: Jacobian of the relative deformation
        //  gradient
        volTensorField relFinv_;

        //- Relative Jacobian: Jacobian of the relative deformation
        //  gradient
        volScalarField relJ_;

        //- Implicit stiffness; coefficient of the Laplacian term
        //  The value of this term only affects convergence and not the
        //  answer In the case of Hookean elasticity, (2*mu + lambda) is
        //  the optimum
        volScalarField impK_;

        //- Implicit stiffness surface field
        surfaceScalarField impKf_;


    // Protected Member Functions

        //- Update fields
        virtual void update(const bool correctSigma = true);

        //- Disallow default bitwise copy construct
        updatedLagSolid(const updatedLagSolid&);

        //- Disallow default bitwise assignment
        void operator=(const updatedLagSolid&);

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::UPDATED_LAGRANGIAN;
        }


public:

    // Constructors

        //- Construct from components
        updatedLagSolid
        (
            const word& type,
            dynamicFvMesh& mesh,
            const bool isSolid = true
        );


    // Destructor

        virtual ~updatedLagSolid()
        {}


    // Member Functions

        // Access

            //- Does this model move the mesh?
            virtual bool movingMesh() const
            {
                return true;
            }


        // Edit
            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;

            //- Update total accumulated fields
            virtual void updateTotalFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "updatedLagSolid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

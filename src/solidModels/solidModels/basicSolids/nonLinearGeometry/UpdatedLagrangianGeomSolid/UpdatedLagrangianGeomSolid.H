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
    UpdatedLagrangianGeomSolid

Description
    Base class for updated Lagrangian solid models

SourceFiles
    UpdatedLagrangianGeomSolid.C

\*---------------------------------------------------------------------------*/

#ifndef UpdatedLagrangianGeomSolid_H
#define UpdatedLagrangianGeomSolid_H

#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                    Class UpdatedLagrangianGeomSolid Declaration
\*---------------------------------------------------------------------------*/

template<class IncrementalModel>
class UpdatedLagrangianGeomSolid
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
        volScalarField relJ_;


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        UpdatedLagrangianGeomSolid(const UpdatedLagrangianGeomSolid&);

        //- Disallow default bitwise assignment
        void operator=(const UpdatedLagrangianGeomSolid&);


public:

    // Constructors

        //- Construct from components
        UpdatedLagrangianGeomSolid
        (
            const word& type,
            dynamicFvMesh& mesh,
            const bool isSolid = true
        );


    // Destructor

        virtual ~UpdatedLagrangianGeomSolid()
        {}


    // Member Functions

        // Access

            //- Return nonlinear geometry enumerator
            virtual nonLinearGeometry::nonLinearType nonLinGeom() const
            {
                return nonLinearGeometry::UPDATED_LAGRANGIAN;
            }

            //- Does this model move the mesh?
            virtual bool movingMesh() const
            {
                return true;
            }

            //- Return face deformation gradient tensor
            virtual tmp<surfaceTensorField> Ff() const
            {
                return fvc::interpolate(F_);
            }

            //- Return Jacobian of a patch
            virtual const scalarField& J(const fvPatch& patch) const
            {
                return J_.boundaryField()[patch.index()];
            }


        // Edit

            //- Update fields
            virtual void update(const bool correctSigma = true);

            //- Return first Piola tensor (relative deformation)
            virtual tmp<volTensorField> P() const;

            //- Return first Piola tensor on a patch
            virtual tmp<tensorField> P(const fvPatch& patch) const;

            //- Boundary normal
            virtual tmp<vectorField> nf(const fvPatch& patch) const;

            //- Update total accumulated fields
            virtual void updateTotalFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "UpdatedLagrangianGeomSolid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

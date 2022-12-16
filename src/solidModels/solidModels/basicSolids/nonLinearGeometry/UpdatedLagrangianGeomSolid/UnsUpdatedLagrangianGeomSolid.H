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
    UnsUpdatedLagrangianGeomSolid

Description
    Base class for unstructured updated Lagrangian solid model

SourceFiles
    UnsUpdatedLagrangianGeomSolid.C

\*---------------------------------------------------------------------------*/

#ifndef UnsUpdatedLagrangianGeomSolid_H
#define UnsUpdatedLagrangianGeomSolid_H

#include "unsSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class UnsUpdatedLagrangianGeomSolid Declaration
\*---------------------------------------------------------------------------*/

template<class IncrementalModel>
class UnsUpdatedLagrangianGeomSolid
:
    public IncrementalModel
{
protected:
    // Protected data

        //- Deformation gradient
        surfaceTensorField Ff_;

        //- Relative deformation gradient
        surfaceTensorField relFf_;

        //- Jacobian of the deformation gradient
        surfaceScalarField Jf_;

        //- Relative Jacobian: Jacobian of the relative deformation gradient
        surfaceScalarField relJf_;

public:

    // Constructors

        //- Construct from components
        UnsUpdatedLagrangianGeomSolid
        (
            const word& type,
            dynamicFvMesh& mesh,
            const bool isSolid = true
        );


    // Destructor

        virtual ~UnsUpdatedLagrangianGeomSolid()
        {}


    // Member Functions

        // Access

            //- Return nonlinear geometry enumerator
            virtual nonLinearGeometry::nonLinearType nonLinGeom() const
            {
                return nonLinearGeometry::UPDATED_LAGRANGIAN;
            }

            //- Return face deformation gradient tensor
            virtual tmp<surfaceTensorField> Ff() const
            {
                return Ff_;
            }

            //- Return Jacobian of a patch
            virtual const scalarField& J(const fvPatch& patch) const
            {
                return Jf_.boundaryField()[patch.index()];
            }


        // Edit

            //- Update fields
            virtual void update(const bool correctSigma = true);

            //- Return first Piola tensor (relative deformation)
            virtual tmp<surfaceTensorField> Pf() const;

            //- Return first Piola tensor on a patch
            virtual tmp<tensorField> Pf(const fvPatch& patch) const;

            //- Boundary normal
            virtual tmp<vectorField> nf(const fvPatch& patch) const;

        //- Write all associated fields
        //  Force write so this should only be used when not
        //  writing at a standard time
        virtual void writeNecessaryFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "UnsUpdatedLagrangianGeomSolid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

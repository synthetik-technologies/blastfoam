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
    TotalLagrangianGeomSolid

Description
    Base class for total Lagrangian solid models

SourceFiles
    TotalLagrangianGeomSolid.C

\*---------------------------------------------------------------------------*/

#ifndef TotalLagrangianGeomSolid_H
#define TotalLagrangianGeomSolid_H

#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class TotalLagrangianGeomSolid Declaration
\*---------------------------------------------------------------------------*/

template<class IncrementalModel>
class TotalLagrangianGeomSolid
:
    public IncrementalModel
{
protected:
    // Protected data

        //- Total deformation gradient
        volTensorField F_;

        //- Relative deformation gradient
        volTensorField relF_;

        //- Jacobian of total deformation gradient
        volScalarField J_;

        //- Relative Jacobian: Jacobian of the relative deformation gradient
        volScalarField relJ_;


public:

    // Constructors

        //- Construct from components
        TotalLagrangianGeomSolid
        (
            const word& type,
            dynamicFvMesh& mesh,
            const bool isSolid = true
        );

    // Destructor

        virtual ~TotalLagrangianGeomSolid()
        {}


    // Member Functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::TOTAL_LAGRANGIAN;
        }

        //- Update fields
        virtual void update(const bool correctSigma = true);

        //- Return first Piola tensor
        virtual tmp<volTensorField> P() const;

        //- Return first Piola tensor on a patch
        virtual tmp<tensorField> P(const fvPatch& patch) const;

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
    #include "TotalLagrangianGeomSolid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

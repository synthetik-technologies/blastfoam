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
    UnsTotalLagrangianGeomSolid

Description
    Base class for unstructured total Lagrangian solid models

SourceFiles
    UnsTotalLagrangianGeomSolid.C

\*---------------------------------------------------------------------------*/

#ifndef UnsTotalLagrangianGeomSolid_H
#define UnsTotalLagrangianGeomSolid_H

#include "unsSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class UnsTotalLagrangianGeomSolid Declaration
\*---------------------------------------------------------------------------*/

template<class IncrementalModel>
class UnsTotalLagrangianGeomSolid
:
    public IncrementalModel
{
protected:
    // Protected data

        //- Total deformation gradient
        surfaceTensorField Ff_;

        //- Inverse of total deformation gradient
        surfaceTensorField invFf_;

        //- Relative deformation gradient
        surfaceTensorField relFf_;

        //- Jacobian of total deformation gradient
        surfaceScalarField Jf_;

        //- Relative Jacobian: Jacobian of the relative deformation gradient
        surfaceScalarField relJf_;


public:

    // Constructors

        //- Construct from components
        UnsTotalLagrangianGeomSolid
        (
            const word& type,
            dynamicFvMesh& mesh,
            const bool isSolid = true
        );

    // Destructor

        virtual ~UnsTotalLagrangianGeomSolid()
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
        virtual tmp<surfaceTensorField> Pf() const;

        //- Return first Piola tensor on a patch
        virtual tmp<tensorField> Pf(const fvPatch& patch) const;

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

        //- Boundary normal
        virtual tmp<vectorField> nf(const fvPatch& patch) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "UnsTotalLagrangianGeomSolid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

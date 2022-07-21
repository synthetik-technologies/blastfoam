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
    unsTotalLagSolid

Description
    Base class for unstructured total Lagrangian solid models

SourceFiles
    unsTotalLagSolid.C

\*---------------------------------------------------------------------------*/

#ifndef unsTotalLagSolid_H
#define unsTotalLagSolid_H

#include "totalLagSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class unsTotalLagSolid Declaration
\*---------------------------------------------------------------------------*/

template<class IncrementalModel>
class unsTotalLagSolid
:
    public totalLagSolid<IncrementalModel>
{
protected:
    // Protected data

        //- Total deformation gradient
        surfaceTensorField Ff_;

        //- Relative deformation gradient
        surfaceTensorField relFf_;

        //- Inverse of the total deformation gradient
        surfaceTensorField Finvf_;

        //- Relative Jacobian: Jacobian of the relative deformation gradient
        surfaceTensorField relFinvf_;

        //- Jacobian of total deformation gradient
        surfaceScalarField Jf_;

        //- Relative Jacobian: Jacobian of the relative deformation gradient
        surfaceScalarField relJf_;


    // Protected member functions

        //- Update fields
        virtual void update(const bool correctSigma = true);


public:

    // Constructors

        //- Construct from components
        unsTotalLagSolid
        (
            const word& type,
            dynamicFvMesh& mesh,
            const bool isSolid = true
        );

    // Destructor

        virtual ~unsTotalLagSolid()
        {}


    // Member Functions

        //- Traction boundary surface normal gradient
        virtual tmp<vectorField> tractionBoundarySnGrad
        (
            const vectorField& traction,
            const scalarField& pressure,
            const fvPatch& patch
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "unsTotalLagSolid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

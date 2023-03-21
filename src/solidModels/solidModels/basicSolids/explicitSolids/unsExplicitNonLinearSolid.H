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
    unsExplicitNonLinearSolid

Description
    Base class for unstructured explicit non linear solid solvers

SourceFiles
    unsExplicitNonLinearSolid.C

\*---------------------------------------------------------------------------*/

#ifndef unsExplicitNonLinearSolid_H
#define unsExplicitNonLinearSolid_H

#include "unsExplicitSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                Class unsExplicitNonLinearSolid Declaration
\*---------------------------------------------------------------------------*/

class unsExplicitNonLinearSolid
:
    public unsExplicitSolid
{
protected:
    // Protected Member Functions

        //- Update wavespeeds
        virtual void updateWavespeeds();

        //- Return first Piola tensor dotted with mesh surface normal
        virtual tmp<vectorField> Pn(const fvPatch& patch) const
        {
            return this->Pf(patch) & patch.nf();
        }

        //- Return face deformation gradient tensor
        virtual tmp<surfaceTensorField> Ff() const = 0;


public:

    // Constructors

        //- Construct from components
        unsExplicitNonLinearSolid
        (
            const word& type,
            dynamicFvMesh& mesh,
            const nonLinearGeometry::nonLinearType nonLinear,
            const bool isSolid
        );

    // Destructor

        virtual ~unsExplicitNonLinearSolid()
        {}


    // Member Functions

        // Access

            //- Return first Piola tensor on a patch
            virtual tmp<tensorField> Pf(const fvPatch& patch) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    LinearGeomSolid

Description
    Base Class for linear solid models

SourceFiles
    LinearGeomSolid.C

\*---------------------------------------------------------------------------*/

#ifndef LinearGeomSolid_H
#define LinearGeomSolid_H

#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class LinearGeomSolid Declaration
\*---------------------------------------------------------------------------*/

template<class IncrementalModel>
class LinearGeomSolid
:
    public IncrementalModel
{
public:

    // Constructors

        //- Construct from components
        LinearGeomSolid
        (
            const word& type,
            fvMesh& mesh,
            const bool isSolid = true
        );

    // Destructor

        virtual ~LinearGeomSolid()
        {}


    // Member Functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }

        //- Update fields
        virtual void update(const bool correctSigma = true);

        //- Boundary normal
        virtual tmp<vectorField> nf(const fvPatch& patch) const
        {
            return patch.nf();
        }

        //- Write
        virtual bool write(const bool write = true) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LinearGeomSolid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

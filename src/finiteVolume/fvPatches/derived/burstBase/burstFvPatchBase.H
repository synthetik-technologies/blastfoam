/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::burstFvPatchBase

Description
    Base class for burst fvPatches

SourceFiles
    burstFvPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef burstFvPatchBase_H
#define burstFvPatchBase_H

#include "cyclicFvPatch.H"
#include "burstPolyPatchBase.H"
#include "fvBoundaryMesh.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class burstFvPatchBase Declaration
\*---------------------------------------------------------------------------*/

class burstFvPatchBase
{
protected:
    // Protected Data

        //- Reference to the fvPatch
        const fvPatch& patch_;

        //- Constant reference to the burstCyclicPolyPatch
        burstPolyPatchBase& burstPolyPatch_;

        //- Constant reference to the intact field
        const scalarField& intact_;

        //- Update the mesh weights and deltas
        void updateDeltas();

        //- Make weights
        virtual void makeWeights(scalarField&) const = 0;


public:

    // Constructors

        //- Construct from polyPatch
        burstFvPatchBase(const fvPatch& patch)
        :
            patch_(patch),
            burstPolyPatch_
            (
                const_cast<burstPolyPatchBase&>
                (
                    refCast<const burstPolyPatchBase>(patch.patch()))
            ),
            intact_(burstPolyPatch_.intact())
        {}

        burstFvPatchBase(const burstFvPatchBase&) = delete;

    // Destructor
    virtual ~burstFvPatchBase()
    {}

    // Member Functions

        //- Return the burstPolyPatchBase
        const burstPolyPatchBase& burstPolyPatch() const
        {
            return burstPolyPatch_;
        }

        //- Access the burstPolyPatchBase
        burstPolyPatchBase& burstPolyPatch()
        {
            return burstPolyPatch_;
        }

        //- Set the polyPatch intact field
        void setIntact(const scalarField& intact);

        //- Update the intact field and update mesh weight and deltas if
        //  a change happens
        void update
        (
            const scalarField& p,
            const scalarField& impulse
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

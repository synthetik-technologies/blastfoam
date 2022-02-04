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
    Foam::burstCyclicFvPatch

Description
    Cyclic-plane patch.

SourceFiles
    burstCyclicFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef burstCyclicFvPatch_H
#define burstCyclicFvPatch_H

#include "cyclicFvPatch.H"
#include "burstCyclicPolyPatch.H"
#include "burstFvPatchBase.H"
#include "fvBoundaryMesh.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class burstCyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

class burstCyclicFvPatch
:
    public cyclicFvPatch,
    public burstFvPatchBase
{
protected:

    // Protected Member functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField& w) const
        {
            cyclicFvPatch::makeWeights(w);
        }


public:

    //- Runtime type information
    TypeName(burstCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        burstCyclicFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            cyclicFvPatch(patch, bm),
            burstFvPatchBase(dynamicCast<const fvPatch>(*this))
        {}


    // Member Functions

        //- Return delta (P to N) vectors across coupled patch
        virtual tmp<vectorField> delta() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

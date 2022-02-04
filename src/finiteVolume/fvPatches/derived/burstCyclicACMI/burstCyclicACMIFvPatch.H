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
    Foam::burstCyclicACMIFvPatch

Description
    Cyclic-plane patch.

SourceFiles
    burstCyclicACMIFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef burstCyclicACMIFvPatch_H
#define burstCyclicACMIFvPatch_H

#include "cyclicACMIFvPatch.H"
#include "burstCyclicACMIPolyPatch.H"
#include "burstFvPatchBase.H"
#include "fvBoundaryMesh.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class burstCyclicACMIFvPatch Declaration
\*---------------------------------------------------------------------------*/

class burstCyclicACMIFvPatch
:
    public cyclicACMIFvPatch,
    public burstFvPatchBase
{
protected:

    // Protected Member functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField& w) const
        {
            cyclicACMIFvPatch::makeWeights(w);
        }

public:

    //- Runtime type information
    TypeName(burstCyclicACMIPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        burstCyclicACMIFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        )
        :
            cyclicACMIFvPatch(patch, bm),
            burstFvPatchBase(dynamicCast<const fvPatch>(*this))
        {}


    // Member Functions

        // Access

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

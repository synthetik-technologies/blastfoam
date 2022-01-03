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
    Foam::burstCyclicAMIFvPatch

Description
    Cyclic-plane patch.

SourceFiles
    burstCyclicAMIFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef burstCyclicAMIFvPatch_H
#define burstCyclicAMIFvPatch_H

#include "cyclicAMIFvPatch.H"
#include "burstCyclicAMIPolyPatch.H"
#include "fvBoundaryMesh.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class burstCyclicAMIFvPatch Declaration
\*---------------------------------------------------------------------------*/

class burstCyclicAMIFvPatch
:
    public cyclicAMIFvPatch
{
    // Private Data

        //- Current time index
        mutable label curTimeIndex_;

        //- Constant reference to the burstCyclicAMIPolyPatch
        const burstCyclicAMIPolyPatch& burstCyclicAMIPolyPatch_;

        //- Intact field
        scalarField& intact_;


public:

    //- Runtime type information
    TypeName(burstCyclicAMIPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        burstCyclicAMIFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            cyclicAMIFvPatch(patch, bm),
            curTimeIndex_(-1),
            burstCyclicAMIPolyPatch_(refCast<const burstCyclicAMIPolyPatch>(patch)),
            intact_
            (
                const_cast<scalarField&>(burstCyclicAMIPolyPatch_.intact())
            )
        {}


    // Member Functions

        // Access

            //- Return local reference cast into the burstCyclicAMI patch
            const burstCyclicAMIPolyPatch& burstCyclicAMIPatch() const
            {
                return burstCyclicAMIPolyPatch_;
            }

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;


            //- Return the intact field
            const Field<scalar>& intact() const
            {
                return intact_;
            }

            //- Access the intact field
            Field<scalar>& intact()
            {
                return intact_;
            }

            bool usePressure() const
            {
                return burstCyclicAMIPolyPatch_.usePressure();
            }

            bool useImpulse() const
            {
                return burstCyclicAMIPolyPatch_.useImpulse();
            }

            //- Return if the patch is up to date
            bool uptoDate(const label curTimeIndex) const
            {
                if (curTimeIndex_ == curTimeIndex)
                {
                    return true;
                }
                curTimeIndex_ = curTimeIndex;
                return false;
            }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatch&, const labelList&);


        // Edit Functions

            //- Update the intact field
            void update(const scalarField& p, const scalarField& impulse);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

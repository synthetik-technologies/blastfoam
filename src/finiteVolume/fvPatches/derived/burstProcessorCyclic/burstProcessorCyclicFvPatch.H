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
    Foam::burstProcessorCyclicFvPatch

Description
    Cyclic-plane patch.

SourceFiles
    burstProcessorCyclicFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef burstProcessorCyclicFvPatch_H
#define burstProcessorCyclicFvPatch_H

#include "cyclicFvPatch.H"
#include "burstProcessorCyclicPolyPatch.H"
#include "fvBoundaryMesh.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class burstProcessorCyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

class burstProcessorCyclicFvPatch
:
    public cyclicFvPatch
{
    // Private Data

        //- Current time index
        mutable label curTimeIndex_;

        //- Constant reference to the burstProcessorCyclicPolyPatch
        const burstProcessorCyclicPolyPatch& burstProcessorCyclicPolyPatch_;

        //- Intact field
        scalarField& intact_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;


public:

    //- Runtime type information
    TypeName(burstProcessorCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        burstProcessorCyclicFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            cyclicFvPatch(patch, bm),
            curTimeIndex_(-1),
            burstProcessorCyclicPolyPatch_(refCast<const burstProcessorCyclicPolyPatch>(patch)),
            intact_
            (
                const_cast<scalarField&>(burstProcessorCyclicPolyPatch_.intact())
            )
        {}


    // Member Functions

        // Access

            //- Return local reference cast into the burstProcessorCyclic patch
            const burstProcessorCyclicPolyPatch& burstProcessorCyclicPatch() const
            {
                return burstProcessorCyclicPolyPatch_;
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
                return burstProcessorCyclicPolyPatch_.usePressure();
            }

            bool useImpulse() const
            {
                return burstProcessorCyclicPolyPatch_.useImpulse();
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

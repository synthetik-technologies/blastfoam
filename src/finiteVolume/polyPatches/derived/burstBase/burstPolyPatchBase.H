/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::burstPolyPatchBase

Description
    Poly patch for a bursting surface (cannot be undone)

SourceFiles
    burstPolyPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef burstPolyPatchBase_H
#define burstPolyPatchBase_H

#include "polyPatch.H"
#include "polyMesh.H"
#include "polyBoundaryMesh.H"
#include "Time.H"
#include "UautoPtr.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                       Class burstPolyPatchBase Declaration
\*---------------------------------------------------------------------------*/


class burstPolyPatchBase
{
protected:
    // Protected Data

        //- Constant reference to the polyPatch
        const polyPatch& patch_;

        //- Intact field
        UautoPtr<const scalarField> intact_;

        //- Is pressure used
        bool usePressure_;

        //- Name of the pressure field
        word pName_;

        //- Burst pressure
        scalar pBurst_;

        //- Reference pressure
        scalar pRef_;

        //- bool impulse used
        bool useImpulse_;

        //- Name of the impulse field
        word impulseName_;

        //- Burst impulse
        scalar impulseBurst_;

        //- Do individual faces burst
        bool partialBurst_;

        //- Intact points
        mutable autoPtr<scalarField> pointIntact_;

        //- Current time index
        mutable label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("burstPolyPatchBase");


    // Constructors

        //- Construct from components
        burstPolyPatchBase
        (
            const polyPatch& p
        );

        //- Construct from dictionary
        burstPolyPatchBase
        (
            const polyPatch& p,
            const dictionary& dict
        );

        //- Construct given the original patch and a map
        burstPolyPatchBase
        (
            const polyPatch& p,
            const burstPolyPatchBase& bppb
        );

        //- Construct given the original patch and a map
        burstPolyPatchBase
        (
            const polyPatch& p,
            const burstPolyPatchBase& bppb,
            const label newSize,
            const label newStart
        );

        //- Construct given the original patch and a map
        burstPolyPatchBase
        (
            const polyPatch& p,
            const burstPolyPatchBase& bppb,
            const labelUList& mapAddressing
        );


    //- Destructor
    virtual ~burstPolyPatchBase();


    // Member Functions

        // Access

            //- Return the polyPatch
            const polyPatch& patch() const
            {
                return patch_;
            }

            //- Access the intact field
            void setIntact(const scalarField& intact)
            {
                if (!intact_.valid())
                {
                    intact_.set(&intact);
                }
            }

            //- Return the interpolated intact field
            const Field<scalar>& pointIntact() const;

            //- Return if pressure is used
            bool usePressure() const
            {
                return usePressure_;
            }

            //- Return the name of the pressure field
            const word& pName() const
            {
                return pName_;
            }

            //- Return burst pressure
            scalar pBurst() const
            {
                return pBurst_;
            }

            //- Return reference pressure
            scalar pRef() const
            {
                return pRef_;
            }

            //- Return if impulse is used
            bool useImpulse() const
            {
                return useImpulse_;
            }

            //- Return the name of the impulse field
            const word& impulseName() const
            {
                return impulseName_;
            }

            //- Return burst impulse
            scalar impulseBurst() const
            {
                return impulseBurst_;
            }

            //- Return if individual faces burst
            bool partialBurst() const
            {
                return partialBurst_;
            }

            //- Return if the patch is up to date
            bool uptoDate(const bool update) const
            {
                label curTimeIndex =
                   patch_.boundaryMesh().mesh().time().timeIndex();

                if (curTimeIndex_ == curTimeIndex)
                {
                    return true;
                }
                if (update)
                {
                    curTimeIndex_ = curTimeIndex;
                }
                return false;
            }

        // Edit Functions

            //- Update the intact field
            bool update
            (
                const scalarField& p,
                const scalarField& impulse,
                scalarField& intact
            ) const;


        // OS

            //- Write the data as a dictionary
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

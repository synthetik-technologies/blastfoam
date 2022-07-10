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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvPatch;
class burstModel;

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

        //- Burst model
        autoPtr<burstModel> burst_;

        //- Intact points
        mutable autoPtr<scalarField> pointIntact_;

        //- Make the pointIntact field
        void makePointIntact() const;


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


        // Edit Functions

            //- Update the intact field
            bool update
            (
                const fvPatch& p,
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

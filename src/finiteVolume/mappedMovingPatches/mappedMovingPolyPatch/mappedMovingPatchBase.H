/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-08-21 Synthetik Applied Technologies: Mapping of patches
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    Foam::mappedMovingPatchBase

Description
    Determines a mapping between patch face centres and processors
    they're on.

SourceFiles
    mappedMovingPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef mappedMovingPatchBase_H
#define mappedMovingPatchBase_H

#include "mappedPatchBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class mappedMovingPatchBase Declaration
\*---------------------------------------------------------------------------*/

class mappedMovingPatchBase
:
    public mappedPatchBase
{

protected:

        //- Displacement field name
        word displacementField_;

        //- Sample patch displacement field name
        word sampleDisplacementField_;

        //- Offset vector (from patch faces to destination mesh objects)
        void setOffsets() const;

        //- Calculate mapping
        virtual void calcMapping() const;

public:

    //- Runtime type information
    TypeName("mappedMovingPatchBase");


    // Constructors

        //- Construct from components
        mappedMovingPatchBase
        (
            const polyPatch& pp
        );

        //- Construct from components
        mappedMovingPatchBase
        (
            const polyPatch& pp,
            const word& sampleRegion,
            const mappedPatchBase::sampleMode mode,
            const word& samplePatch,
            const vectorField& offset
        );

        //- Construct from dictionary
        mappedMovingPatchBase
        (
            const polyPatch& pp,
            const dictionary& dict
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        mappedMovingPatchBase
        (
            const polyPatch& pp,
            const mappedMovingPatchBase& mmpb
        );

        //- Construct given the original patch and a map
        mappedMovingPatchBase
        (
            const polyPatch& pp,
            const mappedMovingPatchBase& mmpb,
            const labelUList& mapAddressing
        );


    //- Destructor
    virtual ~mappedMovingPatchBase();


    // Member Functions

        //- Return the name of the displacement field
        const word& displacementField() const
        {
            return displacementField_;
        }

        //- Return the name of the sample patch displacement field
        const word& sampleDisplacementField() const
        {
            return sampleDisplacementField_;
        }


    // I/O

            //- Write as a dictionary
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

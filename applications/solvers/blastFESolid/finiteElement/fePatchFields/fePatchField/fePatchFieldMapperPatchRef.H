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
    Foam::fePatchFieldMapperPatchRef

Description
    Foam::fePatchFieldMapperPatchRef

\*---------------------------------------------------------------------------*/

#ifndef fePatchFieldMapperPatchRef_H
#define fePatchFieldMapperPatchRef_H

#include "fePatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class fePatchFieldMapperPatchRef Declaration
\*---------------------------------------------------------------------------*/

class fePatchFieldMapperPatchRef
:
    public generalFePatchFieldMapper
{
    // Private Data

        const fePatch& sourcePatch_;
        const fePatch& targetPatch_;


public:

    // Constructors

        //- Construct from components
        fePatchFieldMapperPatchRef
        (
            const fePatch& source,
            const fePatch& target
        )
        :
            sourcePatch_(source),
            targetPatch_(target)
        {}

        //- Disallow default bitwise copy construction
        fePatchFieldMapperPatchRef
        (
            const fePatchFieldMapperPatchRef&
        ) = delete;


    // Member Functions

        const fePatch& sourcePatch() const
        {
            return sourcePatch_;
        }

        const fePatch& targetPatch() const
        {
            return targetPatch_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fePatchFieldMapperPatchRef&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

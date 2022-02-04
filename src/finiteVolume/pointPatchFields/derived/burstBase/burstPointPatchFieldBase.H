/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::burstPointPatchFieldBase

Description

Usage

SourceFiles
    burstPointPatchFieldBase.C

\*---------------------------------------------------------------------------*/

#ifndef burstPointPatchFieldBase_H
#define burstPointPatchFieldBase_H

#include "pointPatchField.H"
#include "pointFields.H"
#include "burstPolyPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class burstPointPatchFieldBase Declaration
\*---------------------------------------------------------------------------*/

class burstPointPatchFieldBase
{
protected:
    // Protected Data

        //- Reference to the burstPolyPatch
        burstPolyPatchBase& burstBase_;


public:

    //- Runtime type information
    TypeName("burstBase");


    // Constructor

        // Construct from polyPatch
        burstPointPatchFieldBase(const polyPatch&);

        //- Construct from polyPatch and a copy
        burstPointPatchFieldBase
        (
            const polyPatch&,
            const burstPointPatchFieldBase& bpf
        );

        //- Construct from polyPatch and mapper
        burstPointPatchFieldBase
        (
            const polyPatch&,
            const burstPointPatchFieldBase& bpf,
            const pointPatchFieldMapper&
        );

    // Destructor
    virtual ~burstPointPatchFieldBase();


    // Member Functions

        //- Return the burst polyPatch base
        const burstPolyPatchBase& burstPolyPatch() const
        {
            return burstBase_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

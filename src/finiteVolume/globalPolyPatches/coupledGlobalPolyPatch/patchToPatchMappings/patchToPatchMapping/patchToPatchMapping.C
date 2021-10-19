/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "patchToPatchMapping.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchToPatchMapping, 0);
    defineRunTimeSelectionTable(patchToPatchMapping, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatchMapping::patchToPatchMapping
(
    const word& type,
    const dictionary& dict,
    const primitivePatch& patchA,
    const primitivePatch& patchB,
    const globalPolyPatch& globalPatchA,
    const globalPolyPatch& globalPatchB
)
:
    dict_(dict.optionalSubDict(type + "Coeffs")),
    patchA_(patchA),
    patchB_(patchB),
    globalPatchA_(globalPatchA),
    globalPatchB_(globalPatchB)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToPatchMapping::~patchToPatchMapping()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::patchToPatchMapping>    Foam::patchToPatchMapping::New
(
    const dictionary& dict,
    const primitivePatch& patchA,
    const primitivePatch& patchB,
    const globalPolyPatch& globalA,
    const globalPolyPatch& globalB
)
{
    const word type(dict.lookup("mappingType"));
    Info<< "Selecting patchToPatchMapping " << type << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown patchToPatchMapping type " << type
            << endl << endl
            << "Valid patchToPatchMapping types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<patchToPatchMapping>
    (
        cstrIter()(dict, patchA, patchB, globalA, globalB)
    );
}


void Foam::patchToPatchMapping::checkFieldSizes
(
    const label fromPatchSize,
    const label toPatchSize,
    const label fromFieldSize,
    const label toFieldSize
) const
{
    if (fromFieldSize != fromPatchSize)
    {
        FatalErrorInFunction
            << "fromField is the wrong size!" << nl
            << "fromField size: " << fromFieldSize
            << ", fromPatch size: " << fromPatchSize
            << abort(FatalError);
    }

    if (toFieldSize != toPatchSize)
    {
        FatalErrorInFunction
            << "toField is wrong size!" << nl
            << "toField size: " << toFieldSize
            << ", toPatch size: " << toPatchSize
            << abort(FatalError);
    }
}


// ************************************************************************* //

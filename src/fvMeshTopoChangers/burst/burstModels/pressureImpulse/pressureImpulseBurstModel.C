/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022-2025
     \\/     M anipulation  | Synthetik Applied Technologies
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

\*---------------------------------------------------------------------------*/

#include "pressureImpulseBurstModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace burstModels
{
    defineTypeNameAndDebug(pressureImpulse, 0);
    addToRunTimeSelectionTable(burstModel, pressureImpulse, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burstModels::pressureImpulse::pressureImpulse
(
    const dictionary& dict,
    const bool coupled
)
:
    field
    (
        dict,
        coupled,
        {
            dict.lookupOrDefault<word>("pName", "p"),
            dict.lookupOrDefault<word>("impulseName", "impulse")
        },
        {
            dict.lookup<scalar>("pBurst"),
            dict.lookup<scalar>("impulseBurst"),
        }
    ),
    pName_(dict.lookupOrDefault<word>("pName", "p")),
    impulseName_(dict.lookupOrDefault<word>("impulseName", "impulse")),
    pRef_(dict.lookupOrDefault<scalar>("pRef", 0.0))
{
    if (!useDelta_)
    {
        burstValues_[pName_] += pRef_;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstModels::pressureImpulse::~pressureImpulse()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

void Foam::burstModels::pressureImpulse::writeData(Ostream& os) const
{
    burstModel::writeData(os);
    writeEntry(os, "pName", pName_);
    writeEntry(os, "pBurst", burstValues_[pName_] - pRef_);
    writeEntry(os, "pRef", pRef_);
    writeEntry(os, "impulseName", impulseName_);
    writeEntry(os, "impulseBurst", burstValues_[impulseName_]);
}


// ************************************************************************* //

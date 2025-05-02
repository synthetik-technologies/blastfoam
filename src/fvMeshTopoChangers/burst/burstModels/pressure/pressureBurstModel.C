/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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

#include "pressureBurstModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace burstModels
{
    defineTypeNameAndDebug(pressure, 0);
    addToRunTimeSelectionTable(burstModel, pressure, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burstModels::pressure::pressure
(
    const dictionary& dict,
    const bool coupled
)
:
    field
    (
        dict,
        coupled,
        {dict.lookupOrDefault<word>("pName", "p")},
        {dict.lookup<scalar>("pBurst")}
    ),
    pRef_(dict.lookupOrDefault<scalar>("pRef", 0.0))
{
    if (!useDelta_)
    {
        burstValues_.begin()() += pRef_;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstModels::pressure::~pressure()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

void Foam::burstModels::pressure::writeData(Ostream& os) const
{
    burstModel::writeData(os);
    writeEntry(os, "pName", burstValues_.begin().key());
    writeEntry(os, "pBurst", burstValues_.begin()() - pRef_);
    writeEntry(os, "pRef", pRef_);
}


// ************************************************************************* //

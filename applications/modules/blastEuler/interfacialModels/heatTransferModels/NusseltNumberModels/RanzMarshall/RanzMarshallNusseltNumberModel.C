/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "RanzMarshallNusseltNumberModel.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace NusseltNumberModels
{
    defineTypeNameAndDebug(RanzMarshall, 0);
    addToRunTimeSelectionTable(NusseltNumberModel, RanzMarshall, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NusseltNumberModels::RanzMarshall::RanzMarshall
(
    const dictionary& dict,
    const phasePair& pair
)
:
    NusseltNumberModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::NusseltNumberModels::RanzMarshall::~RanzMarshall()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::NusseltNumberModels::RanzMarshall::Nu() const
{
    return volScalarField::New
    (
        typeName + "Nu",
        2.0 + 0.6*sqrt(pair_.Re())*cbrt(pair_.Pr())
    );
}


Foam::scalar Foam::NusseltNumberModels::RanzMarshall::cellNu
(
    const label celli
) const
{
    return
        2.0
      + 0.6
       *sqrt(pair_.cellRe(celli))
       *cbrt(pair_.cellPr(celli));
}

// ************************************************************************* //

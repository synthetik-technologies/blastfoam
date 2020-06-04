/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
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

\*---------------------------------------------------------------------------*/

#include "Zuber.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace virtualMassModels
{
    defineTypeNameAndDebug(Zuber, 0);
    addToRunTimeSelectionTable
    (
        virtualMassModel,
        Zuber,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::virtualMassModels::Zuber::Zuber
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    virtualMassModel(dict, pair, registerObject),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict.lookupOrDefault<scalar>
        (
            "residualAlpha",
            pair_.dispersed().residualAlpha().value()
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::virtualMassModels::Zuber::~Zuber()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::virtualMassModels::Zuber::Cvm
(
    const label nodei,
    const label nodej
) const
{

    return
        (1.0 + 2.0*pair_.dispersed().volumeFraction(nodei))
       /max(pair_.continuous().volumeFraction(nodei), residualAlpha_);
}


// ************************************************************************* //

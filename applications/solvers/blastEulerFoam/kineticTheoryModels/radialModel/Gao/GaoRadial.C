/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "GaoRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace radialModels
{
    defineTypeNameAndDebug(Gao, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        Gao,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::Gao::Gao
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    radialModel(dict, kt),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict.lookup("residualAlpha")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::Gao::~Gao()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::Gao::gs0
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    volScalarField g1
    (
        1.0
       /(
           1.0
         - cbrt(min(phase1, kt_.alphaMinFriction())/phase1.alphaMax())
        )
    );
    volScalarField g2
    (
        1.0
       /(
           1.0
         - cbrt(min(phase2, kt_.alphaMinFriction())/phase2.alphaMax())
        )
    );

    return (phase1.d()*g1 + phase2.d()*g2)/(phase1.d() + phase2.d());
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::Gao::gs0prime
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    volScalarField aByaMax
    (
        cbrt(min(max(phase1, scalar(1e-3))/phase1.alphaMax(), 0.999))
    );

    return
        (1.0/(3.0*phase1.alphaMax()))/sqr(aByaMax - sqr(aByaMax))
       *phase1.d()/(phase1.d() + phase2.d());
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "vanderWaals.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace equationOfStates
{
    defineTypeNameAndDebug(vanderWaals, 0);
    addToRunTimeSelectionTable(equationOfState, vanderWaals, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equationOfStates::vanderWaals::vanderWaals
(
    const volScalarField& e,
    const dictionary& dict
)
:
    equationOfState(e, dict),
    a_("a", dimensionSet(-1, 5, -2, 0, 0, 0, 0), dict),
    b_("b", dimensionSet(-1, 3, 0, 0, 0, 0, 0), dict),
    c_
    (
        dimensionedScalar::lookupOrDefault
        (
            "c",
            dict,
            dimensionSet(1, -1, -2, 0, 0, 0, 0),
            0.0
        )
    ),
    gamma_("gamma", dimensionSet(0, 0, 0, 0, 0, 0, 0), dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::equationOfStates::vanderWaals::~vanderWaals()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::equationOfStates::vanderWaals::Gamma() const
{
    return (gamma_ - 1.0)/(1.0 - b_*rho_) + 1.0;
}


Foam::tmp<Foam::volScalarField>
Foam::equationOfStates::vanderWaals::Pi() const
{
    return
        (1.0 - (gamma_ - 1.0)/(1.0 - b_*rho_))*a_*sqr(rho_)
      + ((gamma_ - 1.0)/(1.0 - b_*rho_) + 1.0)*c_;
}


Foam::tmp<Foam::volScalarField>
Foam::equationOfStates::vanderWaals::delta() const
{
    volScalarField p
    (
        (gamma_ - 1.0)/(1.0 - b_*rho_)
       *(rho_*e_ + a_*sqr(rho_) - c_)
      - (a_*sqr(rho_) + c_)
    );
    return
        -b_*(p + a_*sqr(rho_))/(gamma_ - 1.0)
       + ((1.0 - b_*rho_)/(gamma_ - 1.0) - 1.0)*2.0*a_*rho_;
}


// ************************************************************************* //

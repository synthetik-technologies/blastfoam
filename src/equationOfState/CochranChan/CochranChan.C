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

#include "CochranChan.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace equationOfStates
{
    defineTypeNameAndDebug(CochranChan, 0);
    addToRunTimeSelectionTable(equationOfState, CochranChan, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equationOfStates::CochranChan::CochranChan
(
    const volScalarField& e,
    const dictionary& dict
)
:
    equationOfState(e, dict),
    rho0_("rho0", dimensionSet(1, -3, 0, 0, 0, 0, 0), dict),
    e0_("e0", dimensionSet(0, 2, -2, 0, 0, 0, 0), dict),
    Gamma0_("Gamma0", dimless, dict),
    A_("A", dimensionSet(1, -1, -2, 0, 0, 0, 0), dict),
    Epsilon1_("Epsilon1", dimless, dict),
    B_("B", dimensionSet(1, -1, -2, 0, 0, 0, 0), dict),
    Epsilon2_("Epsilon2", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::equationOfStates::CochranChan::~CochranChan()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::equationOfStates::CochranChan::pRef() const
{
    volScalarField rho(stabilizeRho());
    return A_*pow(rho0_/rho, -Epsilon1_) - B_*pow(rho0_/rho, -Epsilon2_);
}


Foam::tmp<Foam::volScalarField>
Foam::equationOfStates::CochranChan::eRef() const
{
    volScalarField rho(stabilizeRho());
    return
      - A_/((1.0 - Epsilon1_)*rho0_)*(pow(rho0_/rho, 1.0 - Epsilon1_) - 1.0)
      + B_/((1.0 - Epsilon2_)*rho0_)*(pow(rho0_/rho, 1.0 - Epsilon2_) - 1.0);
      - e0_;
}


Foam::tmp<Foam::volScalarField>
Foam::equationOfStates::CochranChan::Gamma() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Gamma",
                rho_.time().timeName(),
                rho_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            rho_.mesh(),
            Gamma0_ + 1.0
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::equationOfStates::CochranChan::Pi() const
{
    volScalarField rho(stabilizeRho());
    return Gamma0_*rho_*eRef() - pRef();
}


Foam::tmp<Foam::volScalarField>
Foam::equationOfStates::CochranChan::delta() const
{
    volScalarField rho(stabilizeRho());
    return
        (
           -A_
           *(
               Epsilon1_*pow(rho0_/rho, -Epsilon1_)
              *(Epsilon1_ - Gamma0_ - 1.0)/rho
             + Gamma0_/rho0_
            )/(Epsilon1_ - 1.0)
          + B_
           *(
               Epsilon2_*pow(rho0_/rho, -Epsilon2_)
              *(Epsilon2_ - Gamma0_ - 1.0)/rho
             + Gamma0_/rho0_
            )/(Epsilon2_ - 1.0)
        )/Gamma0_
      - e0_;


}


// ************************************************************************* //

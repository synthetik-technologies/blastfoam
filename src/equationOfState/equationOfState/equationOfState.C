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

#include "equationOfState.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(equationOfState, 0);
    defineRunTimeSelectionTable(equationOfState, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equationOfState::equationOfState
(
    const volScalarField& e,
    const dictionary& dict
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", dict.dictName()),
            e.mesh().time().timeName(),
            e.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        e.mesh(),
        0.0
    ),
    eosDict_(dict),
    e_(e),
    rho_
    (
        IOobject
        (
            IOobject::groupName("rho", dict.dictName()),
            e.mesh().time().timeName(),
            e.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        e.mesh()
    ),
    alphaRho_
    (
        IOobject
        (
            IOobject::groupName("alphaRho", dict.dictName()),
            e.mesh().time().timeName(),
            e.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (*this)*rho_
    ),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", dict.dictName()),
            e.mesh().time().timeName(),
            e.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        e.mesh(),
        dimensionedScalar("0", dimVelocity*dimArea, 0.0)
    ),
    alphaRhoPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi", dict.dictName()),
            e.mesh().time().timeName(),
            e.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        e.mesh(),
        dimensionedScalar("0", dimVelocity*dimArea*dimDensity, 0.0)
    ),
    residualAlpha_("residualAlpha", dimless, dict),
    residualRho_("residualRho", dimDensity, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::equationOfState::~equationOfState()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::equationOfState::e(const volScalarField& p) const
{
    return  (p + Pi())/(Gamma() - 1.0)/stabilizeRho();
}

Foam::tmp<Foam::volScalarField> Foam::equationOfState::pByGamma() const
{
    return rho_*e_ - Pi()/(Gamma() - 1.0);
}


Foam::tmp<Foam::volScalarField> Foam::equationOfState::xi() const
{
    return 1.0/(Gamma() - 1.0);
}


Foam::tmp<Foam::volScalarField>
Foam::equationOfState::cSqr(const volScalarField& P) const
{
    tmp<volScalarField> h((Gamma()*P + Pi())/((Gamma() - 1.0)*stabilizeRho()));
    tmp<volScalarField> xi(1.0/(Gamma() - 1.0));
    return (h - delta())/xi;
}

// ************************************************************************* //

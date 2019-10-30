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

#include "JWL.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace equationOfStates
{
    defineTypeNameAndDebug(JWL, 0);
    addToRunTimeSelectionTable(equationOfState, JWL, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equationOfStates::JWL::JWL
(
    const volScalarField& e,
    const dictionary& dict
)
:
    equationOfState(e, dict),
    active_
    (
        IOobject
        (
            IOobject::groupName("active", rho_.group()),
            rho_.mesh().time().timeName(),
            rho_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho_.mesh(),
        1.0
    ),
    detonationPoints_(),
    vDet_("speed", dimVelocity, great),
    pRef_("pRef", dimPressure, Foam::constant::standard::Pstd.value()),
    rho0_("rho0", dimDensity, dict),
    Gamma0_("Gamma0", dimless, dict),
    A_("A", dimPressure, dict),
    R1_("R1", dimless, dict),
    B_("B", dimPressure, dict),
    R2_("R2", dimless, dict),
    E0_("E0", dimPressure, dict),
    afterburn_(afterburnModel::New(rho_.mesh(), dict))
{
    Switch pointInitiationOnOff = false;
    if (dict.found("initiation"))
    {
        const dictionary& initiationDict(dict.subDict("initiation"));
        pointInitiationOnOff = initiationDict.lookupType<Switch>("active");
        if (pointInitiationOnOff)
        {
            active_ = 0.0;
            detonationPoints_ =
                initiationDict.lookupType<List<vector>>("points");
            vDet_.read(initiationDict);
            pRef_.read(initiationDict);
            Info<< "Initiation Points: " << nl
                << detonationPoints_ << nl
                << "Detonation speed: " << vDet_ << " [m/s]" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::equationOfStates::JWL::~JWL()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::equationOfStates::JWL::update()
{
    if (Foam::min(active_).value() == 1)
    {
        return;
    }

    dimensionedScalar detonationFrontDistance = rho_.time()*vDet_;

    // ----------------------------------------------
    // PV NOTES: add a loop for each initiation point
    // ----------------------------------------------
    forAll(detonationPoints_, pointi)
    {
        dimensionedVector xDet("xDet", dimLength, detonationPoints_[pointi]);
        active_ =
            Foam::max
            (
                active_,
                pos(detonationFrontDistance - mag(rho_.mesh().C() - xDet))
            );
    }

    afterburn_->update();
}


Foam::tmp<Foam::volScalarField>
Foam::equationOfStates::JWL::e(const volScalarField& p) const
{
    if (E0_.value() == 0)
    {
        tmp<volScalarField> eTmp(equationOfState::e(p));

        if (Foam::min(eTmp()*(*this)).value() < 0)
        {
            FatalErrorInFunction
                << "Negative internal energyfor " << group() << " phase." << nl
                << "Either use different pressure and density fields to " << nl
                << "initialize e, or specify a detonation energy, E0." << nl
                << "The detonation energy must be correctly perscribed " << nl
                << "with relation to the JWL coefficients in order to " << nl
                << "obtain physical results. If the detonation energy " << nl
                << "is too low, internal energies can become negative."
                << abort(FatalError);
        };
        return eTmp;
    }

    return E0_/stabilizeRho();
}


Foam::tmp<Foam::volScalarField>
Foam::equationOfStates::JWL::Gamma() const
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
Foam::equationOfStates::JWL::pByGamma() const
{
    return
        (1.0 - active_)*pRef_/Gamma0_
      + active_*equationOfState::pByGamma();
}

Foam::tmp<Foam::volScalarField>
Foam::equationOfStates::JWL::Pi() const
{
    volScalarField V(rho0_/stabilizeRho());
    return
        A_*(Gamma0_/(R1_*V) - 1.0)*exp(-R1_*V)
      + B_*(Gamma0_/(R2_*V) - 1.0)*exp(-R2_*V)
      - Gamma0_*afterburn_->Q()/V;
}


Foam::tmp<Foam::volScalarField>
Foam::equationOfStates::JWL::delta() const
{
    volScalarField rho(stabilizeRho());
    return
        (
            A_*exp(-R1_*rho0_/rho)
           *(Gamma0_*(1.0/(R1_*rho0_) + 1.0/rho) - R1_*rho0_/sqr(rho))
          + B_*exp(-R2_*rho0_/rho)
           *(Gamma0_*(1.0/(R2_*rho0_) + 1.0/rho) - R2_*rho0_/sqr(rho))
        )/Gamma0_
      - afterburn_->Q()/rho0_;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "HuilinPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace granularPressureModels
{
    defineTypeNameAndDebug(Huilin, 0);

    addToRunTimeSelectionTable
    (
        granularPressureModel,
        Huilin,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Huilin::Huilin
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    granularPressureModel(dict, kt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Huilin::~Huilin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Huilin::granularPressure
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& theta1,
    const volScalarField& theta2,
    const volScalarField& g0,
    const dimensionedScalar& e
) const
{
    const scalar pi = Foam::constant::mathematical::pi;
    volScalarField d1(phase1.d());
    volScalarField d2(phase2.d());
    volScalarField d12(0.5*(d1 + d2));

    volScalarField Theta1(theta1);
    Theta1.max(1e-10);
    volScalarField Theta2(theta2);
    Theta2.max(1e-10);

    volScalarField m1(pi/6.0*pow3(d1)*phase1.rho());
    volScalarField m2(pi/6.0*pow3(d2)*phase2.rho());
    volScalarField n1(6.0*phase1/(pi*pow3(d1)));
    volScalarField n2(6.0*phase2/(pi*pow3(d2)));
    volScalarField m0(m1 + m2);
    volScalarField omega
    (
        (m1*Theta1 - m2*Theta2)
       /sqrt
        (
            sqr(m1*Theta1) + sqr(m2*Theta2)
          + Theta1*Theta2*(sqr(m1) + sqr(m2))
        )
    );

    return
        pi*(1 + e)*pow3(d12)*g0*n1*n2*m1*m2*m0*Theta1*Theta2
       /(3.0*(sqr(m1)*Theta1 + sqr(m2)*Theta2))
       *pow
        (
            sqr(m0)*Theta1*Theta2
           /((sqr(m1)*Theta1 + sqr(m2)*Theta2)*(Theta1 + Theta2)),
            3.0/2.0
        )
       *(1.0 - 3.0*omega + 6.0*sqr(omega) - 10.0*pow3(omega));
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Huilin::
granularPressureByAlpha
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& theta1,
    const volScalarField& theta2,
    const volScalarField& g0,
    const volScalarField& g0prime,
    const dimensionedScalar& e
) const
{
    const scalar pi = Foam::constant::mathematical::pi;
    volScalarField d1(phase1.d());
    volScalarField d2(phase2.d());
    volScalarField d12(0.5*(d1 + d2));

    volScalarField Theta1(theta1);
    Theta1.max(1e-10);
    volScalarField Theta2(theta2);
    Theta2.max(1e-10);

    volScalarField m1(pi/6.0*pow3(d1)*phase1.rho());
    volScalarField m2(pi/6.0*pow3(d2)*phase2.rho());
    volScalarField n1(6.0*phase1/(pi*pow3(d1)));
    volScalarField n2(6.0*phase2/(pi*pow3(d2)));
    volScalarField m0(m1 + m2);
    volScalarField omega
    (
        (m1*Theta1 - m2*Theta2)
       /sqrt
        (
            sqr(m1)*sqr(Theta1) + sqr(m2)*sqr(Theta2)
          + Theta1*Theta2*(sqr(m1) + sqr(m2))
        )
    );

    volScalarField granularPressurePrime
    (
        pi*(1.0 + e)*pow3(d12)*n2*m1*m2*m0*Theta1*Theta2
       /(3.0*(sqr(m1)*Theta1 + sqr(m2)*Theta2))
       *pow
        (
            sqr(m0)*Theta1*Theta2
           /((sqr(m1)*Theta1 + sqr(m2)*Theta2)*(Theta1 + Theta2)),
            3.0/2.0
        )
       *(1.0 - 3.0*omega + 6.0*sqr(omega) - 10.0*pow3(omega))
    );

    return 6.0/(pi*pow3(d1))*granularPressurePrime*(g0prime*phase1 + g0);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Huilin::
granularPressureByTheta
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& theta1,
    const volScalarField& theta2,
    const volScalarField& g0,
    const dimensionedScalar& e
) const
{
    const scalar pi = Foam::constant::mathematical::pi;
    volScalarField d1(phase1.d());
    volScalarField d2(phase2.d());
    volScalarField d12(0.5*(d1 + d2));

    volScalarField Theta1(theta1);
    Theta1.max(1e-10);
    volScalarField Theta2(theta2);
    Theta2.max(1e-10);

    volScalarField m1(pi/6.0*pow3(d1)*phase1.rho());
    volScalarField m1Sqr(sqr(m1));
    volScalarField m2(pi/6.0*pow3(d2)*phase2.rho());
    volScalarField m2Sqr(sqr(m2));
    volScalarField n1(6.0*phase1/(pi*pow3(d1)));
    volScalarField n2(6.0*phase2/(pi*pow3(d2)));
    volScalarField m0(m1 + m2);

    volScalarField a(pi*(1.0 + e)*pow3(d12)*g0*n1*n2*m1*m2*m0/3.0);
    tmp<volScalarField> d
    (
        m1Sqr*Theta1*(Theta1 + Theta2) + m2Sqr*Theta2*(Theta1 + Theta2)
    );
    tmp<volScalarField> y
    (
        sqr(m0)*Theta1*Theta2/d()
    );
    volScalarField Y
    (
        pow(y(), 1.5)
    );
    volScalarField delta
    (
        (m1*Theta1 - m2*Theta2)
       /sqrt
        (
            sqr(m1*Theta1)
          + sqr(m2*Theta2)
          + Theta1*Theta2*(m1Sqr + m2Sqr )
        )
    );
    volScalarField Z(1.0 - 3.0*delta + 6.0*sqr(delta) - 10.0*pow3(delta));

    if (&phase1 == &phase2)
    {
        return a/(2.0*sqr(m1))*Y*Z;
    }

    volScalarField X(a*Theta1*Theta2/(m1Sqr*Theta1 + m2Sqr*Theta2));
    tmp<volScalarField> XPrime
    (
        a*sqr(Theta2*m2)/sqr(m1Sqr*Theta1 + m2Sqr*Theta2)
    );
    tmp<volScalarField> dPrime
    (
        2.0*m1Sqr*Theta1 + (m1Sqr + m2Sqr)*Theta2
    );

    tmp<volScalarField> yPrime
    (
        sqr(m0)*Theta2*(1.0/d() - Theta1*dPrime/sqr(d()))
    );
    tmp<volScalarField> YPrime(1.5*yPrime*sqrt(y));

    tmp<volScalarField> n(m1*Theta1 - m2*Theta2);
    tmp<volScalarField> nPrime(m1);
    d.clear();
    d =
        sqr(m1*Theta1) + sqr(m2*Theta2)
      + Theta1*Theta2*(m1Sqr + m2Sqr);
    dPrime.clear();
    dPrime =
        2.0*m1Sqr*Theta1 + Theta2*(m1Sqr + m2Sqr);
    tmp<volScalarField> deltaPrime
    (
        nPrime/sqrt(d()) - 0.5*n*dPrime()/pow(d(), 1.5)
    );
    tmp<volScalarField> ZPrime
    (
        deltaPrime*(3.0 + 12.0*delta - 30.0*sqr(delta))
    );
    return
        XPrime*Y*Z
      + X*YPrime*Z
      + X*Y*ZPrime;
}


// ************************************************************************* //

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
#include "extrapolatedCalculatedFvPatchFields.H"
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
    if (&phase1 == &phase2)
    {
        return volScalarField::New
        (
            "Ps." + phase1.group(),
            sqr(phase1)*phase1.rho()*theta1*2.0*(1.0 + e)*g0
        );
    }

    using Foam::constant::mathematical::pi;

    tmp<volScalarField> tPs
    (
        volScalarField::New
        (
            "Ps." + phase1.group() + "." + phase2.group(),
            phase1.mesh(),
            dimensionedScalar(dimPressure, 0.0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    volScalarField& Ps = tPs.ref();

    tmp<volScalarField> td1(phase1.d());
    const volScalarField& d1 = td1();

    tmp<volScalarField> td2(phase2.d());
    const volScalarField& d2 = td2();

    const volScalarField& rho1 = phase1.rho();
    const volScalarField& rho2 = phase2.rho();

    forAll(Ps, celli)
    {
        const scalar d12 = d1[celli] + d2[celli];
        const scalar m1 = pi/6.0*pow3(d1[celli])*rho1[celli];
        const scalar m2 = pi/6.0*pow3(d2[celli])*rho2[celli];
        const scalar m0 = m1 + m2;
        const scalar n1 = 6.0*phase1[celli]/(pi*pow3(d1[celli]));
        const scalar n2 = 6.0*phase2[celli]/(pi*pow3(d2[celli]));
        const scalar t1 = theta1[celli];
        const scalar t2 = theta2[celli];

        if (t1 > small || t2 > small)
        {
            const scalar omega =
                (m1*t1 - m2*t2)
               /sqrt(sqr(m1*t1) + sqr(m2*t2) + t1*t2*(sqr(m1) + sqr(m2)));

            Ps[celli] =
                pi*(1 + e.value())*pow3(d12)*g0[celli]*n1*n2*m1*m2*m0*t1*t2
               /(3.0*(sqr(m1)*t1 + sqr(m2)*t2))
               *pow
                (
                    sqr(m0)*t1*t2/((sqr(m1)*t1 + sqr(m2)*t2)*(t1 + t2)),
                    1.5
                )
               *(1.0 - 3.0*omega + 6.0*sqr(omega) - 10.0*pow3(omega));
        }
    }
    Ps.correctBoundaryConditions();

    return tPs;
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
    if (&phase1 == &phase2)
    {
        return volScalarField::New
        (
            "dPsdAlpha." + phase1.group(),
            2.0*phase1*phase1.rho()*theta1*(1.0 + e)*(2.0*g0 + phase1*g0prime)
        );
    }

    using Foam::constant::mathematical::pi;

    tmp<volScalarField> td1(phase1.d());
    const volScalarField& d1 = td1();

    tmp<volScalarField> td2(phase2.d());
    const volScalarField& d2 = td2();

    tmp<volScalarField> td12(0.5*(d1 + d2));
    const volScalarField& d12 = td12();

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
    if (&phase1 == &phase2)
    {
        return volScalarField::New
        (
            "dPsdTheta." + phase1.group(),
            sqr(phase1)*phase1.rho()*2.0*(1.0 + e)*g0
        );
    }

    using Foam::constant::mathematical::pi;

    tmp<volScalarField> td1(phase1.d());
    const volScalarField& d1 = td1();

    tmp<volScalarField> td2(phase2.d());
    const volScalarField& d2 = td2();

    tmp<volScalarField> td12(0.5*(d1 + d2));
    const volScalarField& d12 = td12();

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

    volScalarField coeff(pi*(1.0 + e)*pow3(d12)*g0*n1*n2*m1*m2*m0/3.0);
    volScalarField d
    (
        m1Sqr*Theta1*(Theta1 + Theta2) + m2Sqr*Theta2*(Theta1 + Theta2)
    );
    tmp<volScalarField> y
    (
        sqr(m0)*Theta1*Theta2/d
    );
    volScalarField Y
    (
        pow(y(), 1.5)
    );
    volScalarField omega
    (
        (m1*Theta1 - m2*Theta2)
       /sqrt
        (
            sqr(m1*Theta1)
          + sqr(m2*Theta2)
          + Theta1*Theta2*(m1Sqr + m2Sqr )
        )
    );
    volScalarField Z(1.0 - 3.0*omega + 6.0*sqr(omega) - 10.0*pow3(omega));

    volScalarField X(coeff*Theta1*Theta2/(m1Sqr*Theta1 + m2Sqr*Theta2));
    tmp<volScalarField> XPrime
    (
        coeff*sqr(Theta2*m2)/sqr(m1Sqr*Theta1 + m2Sqr*Theta2)
    );
    tmp<volScalarField> dPrime
    (
        2.0*m1Sqr*Theta1 + (m1Sqr + m2Sqr)*Theta2
    );

    tmp<volScalarField> yPrime
    (
        sqr(m0)*Theta2*(1.0/d - Theta1*dPrime/sqr(d))
    );
    tmp<volScalarField> YPrime(1.5*yPrime*sqrt(y));

    tmp<volScalarField> n(m1*Theta1 - m2*Theta2);
    tmp<volScalarField> nPrime(m1);

    d =
        sqr(m1*Theta1) + sqr(m2*Theta2)
      + Theta1*Theta2*(m1Sqr + m2Sqr);
    dPrime.clear();
    dPrime =
        2.0*m1Sqr*Theta1 + Theta2*(m1Sqr + m2Sqr);
    tmp<volScalarField> deltaPrime
    (
        (nPrime*d - 0.5*n*dPrime())/pow(d, 1.5)
    );
    tmp<volScalarField> ZPrime
    (
        deltaPrime*(3.0 + 12.0*omega - 30.0*sqr(omega))
    );
    return
        XPrime*Y*Z
      + X*YPrime*Z
      + X*Y*ZPrime;
}


// ************************************************************************* //

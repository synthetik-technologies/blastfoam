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

#include "HuilinRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace radialModels
{
    defineTypeNameAndDebug(Huilin, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        Huilin,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::Huilin::Huilin
(
    const dictionary& dict,
    const masterSystem& system
)
:
    radialModel(dict, system),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict.lookup("residualAlpha")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::Huilin::~Huilin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::Huilin::gs0
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    tmp<volScalarField> td1 = phase1.d();
    const volScalarField& d1 = td1();

    tmp<volScalarField> td2 = phase2.d();
    const volScalarField& d2 = td2();

    volScalarField f
    (
        "f",
        system_.alphaMax()
       /max
        (
            system_.alphaMax() - system_.alpha(),
            residualAlpha_
        )
    );
    volScalarField delta(0.25*(phase1/d1 + phase2/d2));
    volScalarField dd(d1*d2/(d1 + d2));

    return
        (
            1.0
          + 6.0*dd*delta/f
          + 8.0*sqr(dd*delta*f)
        )*f;
}


Foam::scalar
Foam::kineticTheoryModels::radialModels::Huilin::cellgs0
(
    const label celli,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const scalar d1 = phase1.celld(celli);
    const scalar d2 = phase2.celld(celli);

    const scalar f
    (
        system_.alphaMax()[celli]
       /max
        (
            system_.alphaMax()[celli] - system_.alpha()[celli],
            residualAlpha_
        )
    );
    const scalar delta(0.25*(phase1[celli]/d1 + phase2[celli]/d2));
    const scalar dd(d1*d2/(d1 + d2));

    return
        (
            1.0
          + 6.0*dd*delta*f
          + 8.0*sqr(dd*delta*f)
        )*f;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::Huilin::gs0prime
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    tmp<volScalarField> td1 = phase1.d();
    const volScalarField& d1 = td1();

    tmp<volScalarField> td2 = phase2.d();
    const volScalarField& d2 = td2();

    const volScalarField& alpha = system_.alpha();
    const volScalarField& alphaMax = system_.alphaMax();

    volScalarField f(alphaMax/max(alphaMax() - alpha, residualAlpha_));
    volScalarField fPrime(sqr(f)/alphaMax);
    volScalarField delta(0.25*(phase1/d1 + phase2/d2));
    volScalarField deltaPrime(0.25/d1);
    volScalarField dd(d1*d2/(d1 + d2));

    return
        fPrime
      + 6.0*dd*(deltaPrime*sqr(f) + 2.0*delta*f*fPrime)
      + 8.0*sqr(dd)*(2.0*deltaPrime*delta*pow3(f) + 3.0*sqr(delta*f)*fPrime);
}


Foam::scalar
Foam::kineticTheoryModels::radialModels::Huilin::cellgs0prime
(
    const label celli,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const scalar d1 = phase1.celld(celli);
    const scalar d2 = phase2.celld(celli);

    const scalar alpha = system_.alpha()[celli];
    const scalar alphaMax = system_.alphaMax()[celli];

    const scalar f = alphaMax/max(alphaMax - alpha, residualAlpha_.value());
    const scalar fPrime = sqr(f)/alphaMax;
    const scalar delta(0.25*(phase1[celli]/d1 + phase2[celli]/d2));
    const scalar deltaPrime(0.25/d1);
    const scalar dd(d1*d2/(d1 + d2));

    return
        fPrime
      + 6.0*dd*(deltaPrime*sqr(f) + 2.0*delta*f*fPrime)
      + 8.0*sqr(dd)*(2.0*deltaPrime*delta*pow3(f) + 3.0*sqr(delta*f)*fPrime);
}


// ************************************************************************* //

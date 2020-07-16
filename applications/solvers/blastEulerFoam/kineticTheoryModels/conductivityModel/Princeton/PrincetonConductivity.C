/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "PrincetonConductivity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "phasePairKey.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace conductivityModels
{
    defineTypeNameAndDebug(Princeton, 0);

    addToRunTimeSelectionTable
    (
        conductivityModel,
        Princeton,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::Princeton::Princeton
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    conductivityModel(dict, kt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::Princeton::
~Princeton()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::conductivityModels::Princeton::kappa
(
    const phaseModel& phase,
    const volScalarField& Theta,
    const volScalarField& g0,
    const volScalarField& rho1,
    const volScalarField& da,
    const dimensionedScalar& e
) const
{

    const scalar pi = Foam::constant::mathematical::pi;
    const dimensionedScalar eta = (1.0 + e)/2.0;
    volScalarField alphag0
    (
        IOobject
        (
            "alphag0",
            Theta.time().timeName(),
            Theta.mesh()
        ),
        Theta.mesh(),
        dimensionedScalar("0", dimless, 0.0)
    );

    forAll(kt_.phaseIndexes(), phasei)
    {
        label index2(kt_.phaseIndexes()[phasei]);
        const phaseModel& phase2 = kt_.fluid().phases()[index2];
        phasePairKey key(phase.name(), phase2.name(), false);
//         scalar eij(kt_.es()[key]);
        tmp<volScalarField> gs0ij(kt_.gs0(phase, phase2));

        alphag0 += gs0ij()*phase2;
    }

    tmp<volScalarField> kappa
    (
        75.0*rho1*da*sqrt(pi*Theta)/(48.0*eta*(41.0 - 33.0*eta))
    );
    volScalarField Beta
    (
        IOobject
        (
            "Beta",
            Theta.mesh().time().timeName(),
            Theta.mesh()
        ),
        Theta.mesh(),
        dimensionedScalar("0", dimensionSet(1, -3, -1, 0, 0), 0.0)
    );

    forAllConstIter
    (
        phasePairTable,
        phase.fluid().phasePairs(),
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());
        if
        (
            pair.contains(phase)
         && !kt_.found(pair.otherPhase(phase).name())
        )
        {
            if
            (
                Theta.mesh().foundObject<volScalarField>
                (
                    IOobject::groupName
                    (
                        "Kd",
                        pair.name()
                    )
                )
            )
            {
                Beta += Theta.mesh().template lookupObject<volScalarField>
                (
                    IOobject::groupName
                    (
                        "Kd",
                        pair.name()
                    )
                );
            }
        }
    }

    tmp<volScalarField> kappaStar
    (
        rho1*phase*g0*Theta*kappa()
       /(
            rho1*alphag0*Theta
          + (6.0*Beta*kappa())/(5.0*rho1*max(phase, phase.residualAlpha()))
          + dimensionedScalar("small", dimensionSet(1, -1, -2, 0, 0), SMALL)
        )
    );

    return
    (
        kappaStar/g0
       *(
            (1.0 + 12.0/5.0*eta*alphag0)
           *(1.0 + 12.0/5.0*sqr(eta)*(4.0*eta - 3.0)*alphag0)
          + (
                64.0/(25.0*pi)*(41.0 - 33.0*eta)*sqr(eta)
               *sqr(alphag0)
            )
        )
    );


}


bool Foam::kineticTheoryModels::conductivityModels::Princeton::read()
{
    return true;
}


// ************************************************************************* //

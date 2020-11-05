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

#include "PrincetonViscosity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace viscosityModels
{
    defineTypeNameAndDebug(Princeton, 0);
    addToRunTimeSelectionTable(viscosityModel, Princeton, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::viscosityModels::Princeton::Princeton
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    viscosityModel(dict, kt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::viscosityModels::Princeton::~Princeton()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::viscosityModels::Princeton::nu
(
    const phaseModel& phase,
    const volScalarField& Theta,
    const volScalarField& g0,
    const volScalarField& rho1,
    const volScalarField& da,
    const dimensionedScalar& e
) const
{
    const scalar pi = constant::mathematical::pi;
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

        alphag0 += gs0ij*phase2;
    }

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
    tmp<volScalarField> nu
    (
        5.0/96.0*da*sqrt(pi*Theta)
    );
    tmp<volScalarField> nub
    (
        256.0/(5.0*pi)*nu()*phase*alphag0
    );
    tmp<volScalarField> nuStar
    (
        phase*g0*Theta*nu()
       /(
            max
            (
                alphag0*Theta
              + (2.0*Beta*nu())/(rho1*max(phase, phase.residualAlpha())),
                dimensionedScalar("small", dimensionSet(0, 2, -2, 0, 0), 1e-6)
            )
        )
    );
    return
        3.6/3.0
       *(
            nuStar/(g0*eta*(2.0 - eta))
           *(1.0 + 8.0/5.0*eta*alphag0)
           *(1.0 + 8.0/5.0*eta*(3.0*eta - 2.0)*alphag0)
          + 3.0/5.0*eta*nub
        );
}


// ************************************************************************* //

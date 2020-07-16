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

#include "ChaoConductivity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace conductivityModels
{
    defineTypeNameAndDebug(Chao, 0);

    addToRunTimeSelectionTable
    (
        conductivityModel,
        Chao,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::Chao::Chao
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    conductivityModel(dict, kt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::Chao::
~Chao()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::conductivityModels::Chao::kappa
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

    volScalarField alphaCoeff
    (
        IOobject
        (
            "alphaCoeff",
            Theta.time().timeName(),
            Theta.mesh()
        ),
        Theta.mesh(),
        dimensionedScalar("0", dimless, 0.0)
    );
    volScalarField coeff
    (
        IOobject
        (
            "coeff",
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
        const scalar& eij(kt_.es(phasePairKey(phase.name(), phase2.name())));
        tmp<volScalarField> gs0ij(kt_.gs0(phase, phase2));

        alphaCoeff += gs0ij*phase2*(1.0 + eij);
        coeff += gs0ij*(1.0 + eij);
    }

    coeff /= kt_.phaseIndexes().size();

    return
    (
        150.0/(384.0*coeff)*sqr(1.0 + 6.0/5.0*alphaCoeff)
      + 2.0*phase/pi*alphaCoeff
    )*rho1*da*sqrt(Theta*pi);


}


bool Foam::kineticTheoryModels::conductivityModels::Chao::read()
{
    return true;
}


// ************************************************************************* //

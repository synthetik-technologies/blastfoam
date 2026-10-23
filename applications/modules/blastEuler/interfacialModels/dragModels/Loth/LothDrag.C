/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "LothDrag.H"
#include "phasePair.H"
#include "zeroGradientFvPatchFields.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(Loth, 0);
    addToRunTimeSelectionTable(dragModel, Loth, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Loth::Loth
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dispersedDragModel(dict, pair, registerObject),
    gamma_(dict.lookupOrDefault<scalar>("gamma", -1))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::Loth::~Loth()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::dragModels::Loth::cellGamma(const label celli) const
{
    if (gamma_ > 0)
    {
        return gamma_;
    }

    const blastThermo& thermo = pair_.continuous().thermo();
    const scalar T = thermo.T()[celli];
    return thermo.cellCp(T, celli)/max(thermo.cellCv(T, celli), small);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::dragModels::Loth::CdRe
(
    const scalar Re,
    const scalar Ma,
    const scalar gamma,
    const scalar TpByT
)
{
    using constant::mathematical::pi;

    const scalar Ma4 = pow4(Ma);

    if (Re < 45)
    {
        // Rarefaction dominated
        const scalar Kn = sqrt(0.5*pi*gamma)*Ma/max(Re, small);
        const scalar fKn =
            1.0
           /(
                1.0
              + Kn*(2.514 + 0.8*exp(-0.55/max(Kn, small)))
            );
        const scalar CdReKn = 24.0*(1.0 + 0.15*pow(Re, 0.687))*fKn;

        // Free-molecular term; its weight Ma^4 makes the floor on Ma
        // (which avoids the cancellation in Cd_fm at small s) irrelevant
        const scalar s = max(Ma, 0.01)*sqrt(0.5*gamma);
        const scalar s2 = sqr(s);
        const scalar s4 = sqr(s2);
        const scalar Cdfm =
            (1.0 + 2.0*s2)*exp(-s2)/(s*s2*sqrt(pi))
          + (4.0*s4 + 4.0*s2 - 1.0)*erf(s)/(2.0*s4)
          + 2.0/(3.0*s)*sqrt(pi*max(TpByT, 0.0));
        const scalar CdfmRe =
            Cdfm/(1.0 + (Cdfm/1.63 - 1.0)*sqrt(Re/45.0));

        return (CdReKn + Ma4*CdfmRe*Re)/(1.0 + Ma4);
    }

    // Compression dominated
    const scalar CM =
        Ma < 1.45
      ? 5.0/3.0 + 2.0/3.0*tanh(3.0*log(Ma + 0.1))
      : 2.044 + 0.2*exp(-1.8*sqr(log(Ma/1.5)));
    const scalar GM =
        Ma < 0.89
      ? 1.0 - 1.525*Ma4
      : 0.0002 + 0.0008*tanh(12.77*(Ma - 2.02));
    const scalar HM = 1.0 - 0.258*CM/(1.0 + 514.0*GM);

    return
        24.0*(1.0 + 0.15*pow(Re, 0.687))*HM
      + 0.42*CM*Re/(1.0 + 42500.0*GM/pow(Re, 1.16));
}


Foam::tmp<Foam::volScalarField> Foam::dragModels::Loth::CdRe() const
{
    tmp<volScalarField> tCdRe
    (
        volScalarField::New
        (
            IOobject::groupName("CdRe", pair_.name()),
            pair_.continuous().mesh(),
            dimensionedScalar(dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& cdRe = tCdRe.ref();

    forAll(cdRe, celli)
    {
        cdRe[celli] = cellCdRe(celli);
    }
    cdRe.correctBoundaryConditions();

    return tCdRe;
}


Foam::scalar Foam::dragModels::Loth::cellCdRe(const label celli) const
{
    const phaseModel& continuous = pair_.continuous();
    const scalar p = continuous.cellpRhoT(celli);
    const scalar c = max(continuous.cellSpeedOfSound(p, celli), small);

    return CdRe
    (
        max(pair_.cellRe(celli), small),
        pair_.cellmagUr(celli)/c,
        cellGamma(celli),
        pair_.dispersed().thermo().T()[celli]
       /max(continuous.thermo().T()[celli], small)
    );
}


// ************************************************************************* //

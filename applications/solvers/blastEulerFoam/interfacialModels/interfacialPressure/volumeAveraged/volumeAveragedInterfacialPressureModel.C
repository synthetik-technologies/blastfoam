/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "volumeAveragedInterfacialPressureModel.H"
#include "fluidBlastThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfacialPressureModels
{
    defineTypeNameAndDebug(volumeAveraged, 0);
    addToRunTimeSelectionTable
    (
        interfacialPressureModel,
        volumeAveraged,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacialPressureModels::volumeAveraged::volumeAveraged
(
    const dictionary& dict,
    const phasePair& pair
)
:
    interfacialPressureModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfacialPressureModels::volumeAveraged::~volumeAveraged()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfacialPressureModels::volumeAveraged::PI() const
{
    volScalarField alphaSum(this->pair_.phase1() + this->pair_.phase2());
    alphaSum.max(1e-10);
    volScalarField alphaPI
    (
        this->pair_.phase1()*this->pair_.phase1().p()
      + this->pair_.phase2()*this->pair_.phase2().p()
    );

    return alphaPI/alphaSum;
}


Foam::scalar
Foam::interfacialPressureModels::volumeAveraged::cellPI(const label celli) const
{
    const phaseModel& phase1 = this->pair_.phase1();
    const phaseModel& phase2 = this->pair_.phase2();
    scalar alphaSum(phase1[celli] + phase2[celli]);
    scalar alphaPI
    (
        phase1[celli]*phase1.p()[celli]
      + phase2[celli]*phase2.p()[celli]
    );

    return alphaPI/max(alphaSum, 1e-10);
}


Foam::scalar
Foam::interfacialPressureModels::volumeAveraged::celldPIdAlpha
(
    const label celli,
    const label phasei
) const
{
    const phaseModel& phase1 =
        this->pair_.phase1().index() == phasei
      ? this->pair_.phase1()
      : this->pair_.phase2();
    const phaseModel& phase2 = this->pair_.otherPhase(phase1);
    const scalar p1 = phase1.p()[celli];
    scalar alphaSum
    (
        max(phase1[celli] + phase2[celli], 1e-10)
    );
    scalar alphaPI
    (
        phase1[celli]*p1 + phase2[celli]*phase2.p()[celli]
    );

    return p1/alphaSum - alphaPI/sqr(alphaSum);
}


Foam::scalar
Foam::interfacialPressureModels::volumeAveraged::celldPIde
(
    const label celli,
    const label phasei
) const
{
    const phaseModel& phase1 =
        this->pair_.phase1().index() == phasei
      ? this->pair_.phase1()
      : this->pair_.phase2();
    const phaseModel& phase2 = this->pair_.otherPhase(phase1);
    const fluidBlastThermo& thermo1 =
        dynamicCast<const fluidBlastThermo>(phase1.thermo());
    scalar alphaSum(max(phase1[celli] + phase2[celli], 1e-10));

    return phase1[celli]*thermo1.celldpde(celli)/alphaSum;
}


// ************************************************************************* //

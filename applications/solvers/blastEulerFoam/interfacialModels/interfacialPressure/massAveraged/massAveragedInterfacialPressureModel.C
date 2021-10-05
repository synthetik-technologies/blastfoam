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

#include "massAveragedInterfacialPressureModel.H"
#include "fluidBlastThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfacialPressureModels
{
    defineTypeNameAndDebug(massAveraged, 0);
    addToRunTimeSelectionTable
    (
        interfacialPressureModel,
        massAveraged,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacialPressureModels::massAveraged::massAveraged
(
    const dictionary& dict,
    const phasePair& pair
)
:
    interfacialPressureModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfacialPressureModels::massAveraged::~massAveraged()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfacialPressureModels::massAveraged::PI() const
{
    volScalarField alphaRho
    (
        this->pair_.phase1().alphaRho() + this->pair_.phase2().alphaRho()
    );
    alphaRho.max(1e-10);
    volScalarField alphaRhoPI
    (
        this->pair_.phase1().alphaRho()*this->pair_.phase1().p()
      + this->pair_.phase2().alphaRho()*this->pair_.phase2().p()
    );

    return alphaRhoPI/alphaRho;
}


Foam::scalar
Foam::interfacialPressureModels::massAveraged::cellPI(const label celli) const
{
    const phaseModel& phase1 = this->pair_.phase1();
    const phaseModel& phase2 = this->pair_.phase2();
    scalar alphaRho(phase1.alphaRho()[celli] + phase2.alphaRho()[celli]);
    scalar alphaRhoPI
    (
        phase1.alphaRho()[celli]*phase1.p()[celli]
      + phase2.alphaRho()[celli]*phase2.p()[celli]
    );

    return alphaRhoPI/max(alphaRho, 1e-10);
}


Foam::scalar
Foam::interfacialPressureModels::massAveraged::celldPIdAlpha
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
    const scalar rho1 = phase1.rho()[celli];
    scalar alphaRho
    (
        max(phase1.alphaRho()[celli] + phase2.alphaRho()[celli], 1e-10)
    );
    scalar alphaRhoPI
    (
        phase1.alphaRho()[celli]*p1 + phase2.alphaRho()[celli]*phase2.p()[celli]
    );

    return rho1*p1/alphaRho - alphaRhoPI*rho1/sqr(alphaRho);
}


Foam::scalar
Foam::interfacialPressureModels::massAveraged::celldPIde
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
    scalar alphaRho
    (
        max(phase1.alphaRho()[celli] + phase2.alphaRho()[celli], 1e-10)
    );

    return phase1.alphaRho()[celli]*thermo1.celldpde(celli)/alphaRho;
}


// ************************************************************************* //

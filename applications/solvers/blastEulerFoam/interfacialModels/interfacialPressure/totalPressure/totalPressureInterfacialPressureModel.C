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

#include "totalPressureInterfacialPressureModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfacialPressureModels
{
    defineTypeNameAndDebug(totalPressure, 0);
    addToRunTimeSelectionTable
    (
        interfacialPressureModel,
        totalPressure,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacialPressureModels::totalPressure::totalPressure
(
    const dictionary& dict,
    const phasePair& pair
)
:
    interfacialPressureModel(dict, pair),
    U_(pair.phase1().mesh().lookupObject<volVectorField>("U"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfacialPressureModels::totalPressure::~totalPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfacialPressureModels::totalPressure::PI() const
{
    const phaseModel& phase1 = this->pair_.phase1();
    const phaseModel& phase2 = this->pair_.phase2();
    return
        (
            phase1*(phase1.p() + phase1.rho()*magSqr(phase1.U() - U_))
          + phase2*(phase2.p() + phase2.rho()*magSqr(phase2.U() - U_))
        )/max(phase1 + phase2, 1e-10);
}


Foam::scalar
Foam::interfacialPressureModels::totalPressure::cellPI(const label celli) const
{
    const phaseModel& phase1 = this->pair_.phase1();
    const phaseModel& phase2 = this->pair_.phase2();
    return
        (
            phase1[celli]
           *(
                phase1.p()[celli]
              + phase1.rho()[celli]*magSqr(phase1.U()[celli] - U_[celli])
            )
          + phase2[celli]
           *(
                phase2.p()[celli]
              + phase2.rho()[celli]*magSqr(phase2.U()[celli] - U_[celli])
            )
        )/max(phase1[celli] + phase2[celli], 1e-10);
}


Foam::scalar
Foam::interfacialPressureModels::totalPressure::celldPIdAlpha
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
    const scalar K1 = magSqr(phase1.U()[celli] - U_[celli]);
    const scalar sumAlpha(max(phase1[celli] + phase2[celli], 1e-10));
    scalar pI
    (
        (
            phase1[celli]*(p1 + rho1*K1)
          + phase2[celli]
           *(
                phase2.p()[celli]
              + phase2.rho()[celli]*magSqr(phase2.U()[celli] - U_[celli])
            )
        )/sumAlpha
    );

    return p1/sumAlpha - pI/sqr(sumAlpha);
}


Foam::scalar
Foam::interfacialPressureModels::totalPressure::celldPIde
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
    scalar sumAlpha
    (
        max(phase1[celli] + phase2[celli], 1e-10)
    );

    return phase1[celli]*thermo1.celldpde(celli)/sumAlpha;
}


// ************************************************************************* //

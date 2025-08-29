/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "noPhaseChange.H"
#include "fvScalarMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeModels
{
    defineTypeNameAndDebug(noPhaseChange, 0);
    addToRunTimeSelectionTable(phaseChangeModel, noPhaseChange, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModels::noPhaseChange::noPhaseChange
(
    const dictionary& dict,
    const volScalarField& alpha1,
    const phaseFluidBlastThermo& thermo1,
    const volScalarField& alpha2,
    const phaseFluidBlastThermo& thermo2
)
:
    phaseChangeModel(typeName, dict, alpha1, thermo1, alpha2, thermo2)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::phaseChangeModels::noPhaseChange::mDots() const
{
    return Pair<tmp<volScalarField>>
    (
        volScalarField::null(),
        volScalarField::null()
    );
}


void Foam::phaseChangeModels::noPhaseChange::correct()
{
    phaseChangeModel::correct();
}


// ************************************************************************* //

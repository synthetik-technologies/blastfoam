/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "Kunz.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeModels
{
    defineTypeNameAndDebug(Kunz, 0);
    addToRunTimeSelectionTable(cavitationModel, Kunz, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModels::Kunz::Kunz
(
    const dictionary& dict,
    const volScalarField& alpha1,
    const phaseFluidBlastThermo& thermo1,
    const volScalarField& alpha2,
    const phaseFluidBlastThermo& thermo2
)
:
    cavitationModel(typeName, dict, alpha1, thermo1, alpha2, thermo2),

    UInf_("UInf", dimVelocity, dict_),
    tInf_("tInf", dimTime, dict_),
    Cc_("Cc", dimless, dict_),
    Cv_("Cv", dimless, dict_),

    mcCoeff_(Cc_/tInf_),
    mvCoeff_(Cv_/(0.5*sqr(UInf_)*tInf_)),

    minRhol_("minRhol", dimDensity, dict)
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::phaseChangeModels::Kunz::mDotCV() const
{
    const volScalarField::Internal& p =
        alpha1_.db().lookupObject<volScalarField>("p");

    return Pair<tmp<volScalarField::Internal>>
    (
        mcCoeff_*alphav()*rhov()*sqr(alphal())*pos0(p - pSat()),
        mvCoeff_*alphal()*rhov()/max(rhol(), minRhol_)*max(pSat() - p, p0())
    );
}


void Foam::phaseChangeModels::Kunz::correct()
{
    cavitationModel::correct();
}


// ************************************************************************* //

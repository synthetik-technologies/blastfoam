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

#include "cavitationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeModels
{
    defineTypeNameAndDebug(cavitationModel, 0);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModels::cavitationModel::cavitationModel
(
    const word& type,
    const dictionary& dict,
    const volScalarField& alpha1,
    const phaseFluidBlastThermo& thermo1,
    const volScalarField& alpha2,
    const phaseFluidBlastThermo& thermo2
)
:
    phaseChangeModel(type, dict, alpha1, thermo1, alpha2, thermo2),
    phase1Liquid_(false),
    pSat_("pSat", dimPressure, dict.lookup("pSat")),
    p0_("0", pSat().dimensions(), 0.0)
{
    const word liquidPhaseName(dict.lookup<word>("liquid"));
    if (liquidPhaseName == alpha1.group())
    {
        phase1Liquid_ = true;
    }
    else if (liquidPhaseName != alpha2.group())
    {
        FatalIOErrorInFunction(dict)
            << "Liquid phase is neither " << alpha1.group() << " or "
            << alpha2.group() << endl
            << abort(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::phaseChangeModels::cavitationModel::mDots() const
{
    Pair<tmp<volScalarField::Internal>> dmdts(this->mDotCV());
    if (!phase1Liquid_)
    {
        return reverse(dmdts);
    }
    return dmdts;
}

// ************************************************************************* //

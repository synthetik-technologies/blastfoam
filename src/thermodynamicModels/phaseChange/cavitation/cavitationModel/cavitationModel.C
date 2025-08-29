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
    defineTypeNameAndDebug(cavitationModel, 0);
    defineRunTimeSelectionTable(cavitationModel, phaseThermo);
    defineRunTimeSelectionTable(cavitationModel, thermo);

    const dimensionedScalar cavitationModel::zeroP(dimPressure, 0.0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cavitationModel::cavitationModel
(
    const word& type,
    const dictionary& dict,
    const volScalarField& alpha1,
    const volScalarField& rho1,
    const volScalarField& T1,
    const volScalarField& alpha2,
    const volScalarField& rho2,
    const volScalarField& T2
)
:
    phaseChangeModel(type, dict, alpha1, alpha2),
    phase1Liquid_(false),
    rho1_(rho1),
    rho2_(rho2),
    T1_(T1),
    T2_(T2),
    pSat_(saturationPressureModel::New("pSat", dict_))
{
    const word liquidPhaseName(dict_.lookup<word>("liquid"));
    if (liquidPhaseName == alpha1.group())
    {
        phase1Liquid_ = true;
    }
    else if (liquidPhaseName != alpha2.group())
    {
        FatalIOErrorInFunction(dict_)
            << "Liquid phase is neither " << alpha1.group() << " or "
            << alpha2.group() << endl
            << abort(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::cavitationModel::mDots() const
{
    Pair<tmp<volScalarField::Internal>> dmdts(this->mDotCV());
    if (!phase1Liquid_)
    {
        return reverse(dmdts);
    }
    return dmdts;
}

// ************************************************************************* //

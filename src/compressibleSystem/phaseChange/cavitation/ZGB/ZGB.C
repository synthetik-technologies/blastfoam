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

#include "ZGB.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeModels
{
    defineTypeNameAndDebug(ZGB, 0);
    addToRunTimeSelectionTable
    (
        cavitationModel,
        ZGB,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModels::ZGB::ZGB
(
    const dictionary& dict,
    const volScalarField& alpha1,
    const phaseFluidBlastThermo& thermo1,
    const volScalarField& alpha2,
    const phaseFluidBlastThermo& thermo2
)
:
    cavitationModel(typeName, dict, alpha1, thermo1, alpha2, thermo2),

    n_("n", dimless/dimVolume, dict_),
    alphaNuc_("alphaNuc", dimless, 0.0),
    Cc_("Cc", dimless, dict_),
    Cv_("Cv", dimless, dict_),

    residualRho_
    (
        "residualRho",
        max(thermo1_.residualRho(), thermo2_.residualRho())
    ),
    minRhol_("minRhol", dimDensity, dict)
{
    if (dict_.found("alphaNuc"))
    {
        alphaNuc_ = dict.lookup<scalar>("alphaNuc");
    }
    else
    {
        dimensionedScalar d("dNuc", dimLength, dict_);
        alphaNuc_ =
            (n_*constant::mathematical::pi*pow3(d)/6.0)
           /(1.0 + n_*constant::mathematical::pi*pow3(d)/6.0);
    }
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::phaseChangeModels::ZGB::rRb
(
    const volScalarField::Internal& alphal
) const
{
    using Foam::constant::mathematical::pi;
    return cbrt((4.0/3.0)*pi*n_*alphal/(1.0 + alphaNuc_ - alphal));
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::phaseChangeModels::ZGB::mDotCV() const
{
    const volScalarField::Internal& p =
        alpha1_.db().lookupObject<volScalarField>("p");
    const volScalarField::Internal& alphav = this->alphav();
    const volScalarField::Internal& alphal = this->alphal();
    const volScalarField::Internal& rhov = this->rhov();
    const volScalarField::Internal& rhol = this->rhol();

    const volScalarField::Internal pCoeff
    (
        "pCoeff",
        3.0*rhov*rRb(alphal)
       *sqrt(2.0/(3.0*max(rhol, minRhol_)))
    );

    return Pair<tmp<volScalarField::Internal>>
    (
        Cc_*alphac*pCoeff*sqrt(max(p - pSat(), p0())),
        Cv_*alphal*alphaNuc_*pCoeff*sqrt(max(pSat() - p, p0()))
    );
}


void Foam::phaseChangeModels::ZGB::correct()
{
    cavitationModel::correct();
}


// ************************************************************************* //

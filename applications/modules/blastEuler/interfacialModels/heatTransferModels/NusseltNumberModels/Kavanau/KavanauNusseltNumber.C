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

#include "KavanauNusseltNumber.H"
#include "phasePair.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace NusseltNumberModels
{
    defineTypeNameAndDebug(Kavanau, 0);
    addToRunTimeSelectionTable(NusseltNumberModel, Kavanau, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NusseltNumberModels::Kavanau::Kavanau
(
    const dictionary& dict,
    const phasePair& pair
)
:
    NusseltNumberModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::NusseltNumberModels::Kavanau::~Kavanau()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::NusseltNumberModels::Kavanau::Nu
(
    const scalar Re,
    const scalar Pr,
    const scalar MaByRe
)
{
    const scalar Nuc = 2.0 + 0.459*pow(max(Re, 0.0), 0.55)*pow(Pr, 0.33);
    return Nuc/(1.0 + 3.42*MaByRe*Nuc/max(Pr, small));
}


Foam::tmp<Foam::volScalarField>
Foam::NusseltNumberModels::Kavanau::Nu() const
{
    tmp<volScalarField> tNu
    (
        volScalarField::New
        (
            IOobject::groupName("Nu", pair_.name()),
            pair_.continuous().mesh(),
            dimensionedScalar(dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& nu = tNu.ref();

    forAll(nu, celli)
    {
        nu[celli] = cellNu(celli);
    }
    nu.correctBoundaryConditions();

    return tNu;
}


Foam::scalar Foam::NusseltNumberModels::Kavanau::cellNu
(
    const label celli
) const
{
    const phaseModel& continuous = pair_.continuous();
    const scalar c =
        max
        (
            continuous.cellSpeedOfSound(continuous.cellpRhoT(celli), celli),
            small
        );
    const scalar d = max(pair_.dispersed().celld(celli), 1e-10);

    return Nu
    (
        pair_.cellRe(celli),
        pair_.cellPr(celli),
        max(continuous.cellnu(celli), 1e-15)/(c*d)
    );
}


// ************************************************************************* //

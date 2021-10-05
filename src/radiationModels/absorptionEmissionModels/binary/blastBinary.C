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

#include "blastBinary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(blastBinary, 0);

    addToRunTimeSelectionTable
    (
        blastAbsorptionEmissionModel,
        blastBinary,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::blastBinary::blastBinary
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    blastAbsorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs")),
    model1_
    (
        blastAbsorptionEmissionModel::New(coeffsDict_.subDict("model1"), mesh)
    ),
    model2_
    (
        blastAbsorptionEmissionModel::New(coeffsDict_.subDict("model2"), mesh)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::blastBinary::~blastBinary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::blastBinary::aCont
(
    const label bandI
) const
{
    return model1_->aCont(bandI) + model2_->aCont(bandI);
}


Foam::scalar
Foam::radiationModels::absorptionEmissionModels::blastBinary::cellaCont
(
    const label celli,
    const label bandI
) const
{
    return
        model1_->cellaCont(celli, bandI)
      + model2_->cellaCont(celli, bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::blastBinary::aDisp
(
    const label bandI
) const
{
    return model1_->aDisp(bandI) + model2_->aDisp(bandI);
}



Foam::scalar
Foam::radiationModels::absorptionEmissionModels::blastBinary::cellaDisp
(
    const label celli,
    const label bandI
) const
{
    return
        model1_->cellaDisp(celli, bandI)
      + model2_->cellaDisp(celli, bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::blastBinary::eCont
(
    const label bandI
) const
{
    return model1_->eCont(bandI) + model2_->eCont(bandI);
}


Foam::scalar
Foam::radiationModels::absorptionEmissionModels::blastBinary::celleCont
(
    const label celli,
    const label bandI
) const
{
    return
        model1_->celleCont(celli, bandI)
      + model2_->celleCont(celli, bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::blastBinary::eDisp
(
    const label bandI
) const
{
    return model1_->eDisp(bandI) + model2_->eDisp(bandI);
}


Foam::scalar
Foam::radiationModels::absorptionEmissionModels::blastBinary::celleDisp
(
    const label celli,
    const label bandI
) const
{
    return
        model1_->celleDisp(celli, bandI)
      + model2_->celleDisp(celli, bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::blastBinary::ECont
(
    const label bandI
) const
{
    return model1_->ECont(bandI) + model2_->ECont(bandI);
}


Foam::scalar
Foam::radiationModels::absorptionEmissionModels::blastBinary::cellECont
(
    const label celli,
    const label bandI
) const
{
    return
        model1_->cellECont(celli, bandI)
      + model2_->celleCont(celli, bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::blastBinary::EDisp
(
    const label bandI
) const
{
    return model1_->EDisp(bandI) + model2_->EDisp(bandI);
}


Foam::scalar
Foam::radiationModels::absorptionEmissionModels::blastBinary::cellEDisp
(
    const label celli,
    const label bandI
) const
{
    return
        model1_->cellEDisp(celli, bandI)
      + model2_->cellEDisp(celli, bandI);
}

// ************************************************************************* //

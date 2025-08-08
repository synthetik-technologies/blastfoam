/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "noCombustionBlastCombustionModel.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(noCombustion, 0);
    addToRunTimeSelectionTable(blastCombustionModel, noCombustion, fluid);
    addToRunTimeSelectionTable(blastCombustionModel, noCombustion, solid);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::noCombustion::noCombustion
(
    const word& modelType,
    const multicomponentBlastThermo& thermo,
    const bool isFluid,
    const word& combustionProperties
)
:
    blastCombustionModel(modelType, thermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::noCombustion::~noCombustion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::noCombustion::correct()
{}


Foam::tmp<Foam::volScalarField::Internal>
Foam::combustionModels::noCombustion::R(const label speciei) const
{
    return
        volScalarField::Internal::New
        (
            typedName("R_" + this->thermo().Y()[speciei].name()),
            this->mesh(),
            dimensionedScalar(dimDensity/dimTime, 0)
        );
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::noCombustion::R(volScalarField& Y) const
{
    return tmp<fvScalarMatrix>(new fvScalarMatrix(Y, dimMass/dimTime));
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::noCombustion::Qdot() const
{
    return volScalarField::New
    (
        this->thermo().phasePropertyName(typedName("Qdot")),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    );
}


bool Foam::combustionModels::noCombustion::read()
{
    if (blastCombustionModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technologies
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "constantMassDiameterModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(constantMassDiameter, 0);
    addToRunTimeSelectionTable(diameterModel, constantMassDiameter, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::constantMassDiameter::constantMassDiameter
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    diameterModel(mesh, dict, phaseName),
    rho0_("rho0", dimDensity, dict),
    d0_("d0", dimLength, dict),
    M0_(Foam::constant::mathematical::pi/6.0*pow3(d0_)*rho0_)
{
    const volScalarField& rho =
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("rho", phaseName)
        );
    this->d_ = cbrt(M0_*6.0/(max(rho, 1e-10)*Foam::constant::mathematical::pi));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::constantMassDiameter::~constantMassDiameter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::constantMassDiameter::solve()
{
    const volScalarField& rho =
        this->d_.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("rho", this->d_.group())
        );
    this->d_ = cbrt(M0_*6.0/(max(rho, 1e-10)*Foam::constant::mathematical::pi));
}

// ************************************************************************* //

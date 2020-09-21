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

#include "massAveragedInterfacialVelocityModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfacialVelocityModels
{
    defineTypeNameAndDebug(massAveraged, 0);
    addToRunTimeSelectionTable
    (
        interfacialVelocityModel,
        massAveraged,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacialVelocityModels::massAveraged::massAveraged
(
    const dictionary& dict,
    const phasePair& pair
)
:
    interfacialVelocityModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfacialVelocityModels::massAveraged::~massAveraged()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::interfacialVelocityModels::massAveraged::Ui() const
{
    volScalarField alphaRho
    (
        this->pair_.phase1().alphaRho() + this->pair_.phase2().alphaRho()
    );
    volVectorField alphaRhoUi
    (
        this->pair_.phase1().alphaRho()*this->pair_.phase1().U()
      + this->pair_.phase2().alphaRho()*this->pair_.phase2().U()
    );

    return alphaRhoUi/max(alphaRho, dimensionedScalar(dimDensity, 1e-10));
}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfacialVelocityModels::massAveraged::phi() const
{
    volScalarField alphaRho
    (
        this->pair_.phase1().alphaRho() + this->pair_.phase2().alphaRho()
    );
    surfaceScalarField alphaRhoPhii
    (
        this->pair_.phase1().alphaRhoPhi()
      + this->pair_.phase2().alphaRhoPhi()
    );

    return
        alphaRhoPhii
       /fvc::interpolate(max(alphaRho, dimensionedScalar(dimDensity, 1e-10)));
}


// ************************************************************************* //

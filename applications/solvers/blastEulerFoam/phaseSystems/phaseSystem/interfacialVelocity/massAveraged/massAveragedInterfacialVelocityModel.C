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
    const phaseModelList& phaseModels
)
:
    interfacialVelocityModel(dict, phaseModels),
    phaseNames_(dict.lookup("phases"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfacialVelocityModels::massAveraged::~massAveraged()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::interfacialVelocityModels::massAveraged::Ui() const
{
    volScalarField alphaRho(phaseModels_[phaseNames_[0]].alphaRho());
    volVectorField alphaRhoUi(phaseModels_[phaseNames_[0]].alphaRhoU());
    for (label i = 1; i < phaseNames_.size(); i++)
    {
        alphaRho += phaseModels_[phaseNames_[i]].alphaRho();
        alphaRhoUi += (phaseModels_[phaseNames_[i]].alphaRhoU());
    }

    return alphaRhoUi/max(alphaRho, dimensionedScalar(dimDensity, 1e-10));
}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfacialVelocityModels::massAveraged::phi() const
{
    surfaceScalarField alphaRhoPhi
    (
        phaseModels_[phaseNames_[0]].alphaRhoPhi()
    );
    volScalarField alphaRho(phaseModels_[phaseNames_[0]].alphaRho());
    for (label i = 1; i < phaseNames_.size(); i++)
    {
        alphaRhoPhi += phaseModels_[phaseNames_[i]].alphaRhoPhi();
        alphaRho += (phaseModels_[phaseNames_[i]].alphaRho());
    }

    return
        alphaRhoPhi
       /max
        (
            fvc::interpolate(alphaRho),
            dimensionedScalar(dimDensity, 1e-10)
        );
}


// ************************************************************************* //

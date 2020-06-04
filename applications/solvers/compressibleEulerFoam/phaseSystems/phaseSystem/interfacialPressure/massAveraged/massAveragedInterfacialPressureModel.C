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

#include "massAveragedInterfacialPressureModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfacialPressureModels
{
    defineTypeNameAndDebug(massAveraged, 0);
    addToRunTimeSelectionTable
    (
        interfacialPressureModel,
        massAveraged,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacialPressureModels::massAveraged::massAveraged
(
    const dictionary& dict,
    const phaseModelList& phaseModels
)
:
    interfacialPressureModel(dict, phaseModels),
    phaseNames_(dict.lookup("phases"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfacialPressureModels::massAveraged::~massAveraged()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfacialPressureModels::massAveraged::Pi() const
{
    volScalarField alphaRho(phaseModels_[phaseNames_[0]].alphaRho());
    volScalarField alphaRhoPi(alphaRho*phaseModels_[phaseNames_[0]].p());
    for (label i = 1; i < phaseNames_.size(); i++)
    {
        const volScalarField& alphaRhoi =
            phaseModels_[phaseNames_[i]].alphaRho();
        alphaRho += alphaRhoi;
        alphaRhoPi += (alphaRhoi*phaseModels_[phaseNames_[i]].p());
    }

    return alphaRhoPi/max(alphaRho, dimensionedScalar(dimDensity, 1e-10));;
}


// ************************************************************************* //

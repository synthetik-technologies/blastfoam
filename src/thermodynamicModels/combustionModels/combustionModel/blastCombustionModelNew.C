/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "blastCombustionModel.H"
#include "noCombustionBlastCombustionModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::blastCombustionModel>
Foam::blastCombustionModel::NewFluid
(
    const multicomponentBlastThermo& thermo,
    const word& combustionProperties,
    bool required
)
{
    typeIOobject<IOdictionary> combIO
    (
        IOobject
        (
            thermo.phasePropertyName(combustionProperties),
            thermo.T().mesh().time().constant(),
            thermo.T().mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word modelType(combustionModels::noCombustion::typeName);
    if (combIO.headerOk() || required)
    {
        IOdictionary(combIO).lookup(blastCombustionModel::typeName) >> modelType;
    }
    else
    {
        Info<< "Combustion model not active: "
            << thermo.phasePropertyName(combustionProperties)
            << " not found" << endl;
    }

    Info<< "Selecting combustion model " << modelType << endl;

    const wordList cmpts2(basicThermo::splitThermoName(modelType, 2));
    const wordList cmpts3(basicThermo::splitThermoName(modelType, 3));
    if (cmpts2.size() == 2 || cmpts3.size() == 3)
    {
        modelType = cmpts2.size() ? cmpts2[0] : cmpts3[0];

        WarningInFunction
            << "Template parameters are no longer required when selecting a "
            << blastCombustionModel::typeName << ". This information is now "
            << "obtained directly from the thermodynamics. Actually selecting "
            << "combustion model " << modelType << "." << endl;
    }

    typename fluidConstructorTable::iterator cstrIter =
        fluidConstructorTablePtr_->find(modelType);

    if (cstrIter == fluidConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fluid " << blastCombustionModel::typeName << " type "
            << modelType << nl << nl
            << "Valid " << blastCombustionModel::typeName << " types are:" << nl
            << fluidConstructorTablePtr_->sortedToc()
            << exit(FatalError);

        const wordList names(fluidConstructorTablePtr_->sortedToc());
    }

    return autoPtr<blastCombustionModel>
    (
        cstrIter()(modelType, thermo, true, combustionProperties)
    );
}


Foam::autoPtr<Foam::blastCombustionModel>
Foam::blastCombustionModel::NewSolid
(
    const multicomponentBlastThermo& thermo,
    const word& combustionProperties,
    bool required
)
{
    typeIOobject<IOdictionary> combIO
    (
        IOobject
        (
            thermo.phasePropertyName(combustionProperties),
            thermo.T().mesh().time().constant(),
            thermo.T().mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word modelType(combustionModels::noCombustion::typeName);
    if (combIO.headerOk() || required)
    {
        IOdictionary(combIO).lookup(blastCombustionModel::typeName) >> modelType;
    }
    else
    {
        Info<< "Combustion model not active: "
            << thermo.phasePropertyName(combustionProperties)
            << " not found" << endl;
    }

    Info<< "Selecting combustion model " << modelType << endl;

    const wordList cmpts2(basicThermo::splitThermoName(modelType, 2));
    const wordList cmpts3(basicThermo::splitThermoName(modelType, 3));
    if (cmpts2.size() == 2 || cmpts3.size() == 3)
    {
        modelType = cmpts2.size() ? cmpts2[0] : cmpts3[0];

        WarningInFunction
            << "Template parameters are no longer required when selecting a "
            << blastCombustionModel::typeName << ". This information is now "
            << "obtained directly from the thermodynamics. Actually selecting "
            << "combustion model " << modelType << "." << endl;
    }

    typename solidConstructorTable::iterator cstrIter =
        solidConstructorTablePtr_->find(modelType);

    if (cstrIter == solidConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown solid " << blastCombustionModel::typeName << " type "
            << modelType << nl << nl
            << "Valid " << blastCombustionModel::typeName << " types are:" << nl
            << solidConstructorTablePtr_->sortedToc()
            << exit(FatalError);

        const wordList names(solidConstructorTablePtr_->sortedToc());
    }

    return autoPtr<blastCombustionModel>
    (
        cstrIter()(modelType, thermo, false, combustionProperties)
    );
}


// ************************************************************************* //

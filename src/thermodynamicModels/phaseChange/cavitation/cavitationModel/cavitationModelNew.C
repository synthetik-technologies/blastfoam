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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::cavitationModel> Foam::cavitationModel::New
(
    const dictionary& dict,
    const volScalarField& alpha1,
    const phaseFluidBlastThermo& thermo1,
    const volScalarField& alpha2,
    const phaseFluidBlastThermo& thermo2
)
{
    word modelType(dict.lookup(typeName));

    Info<< "Selecting " << typeName << ": " << modelType << endl;

    phaseThermoConstructorTable::iterator cstrIter =
        phaseThermoConstructorTablePtr_->find(modelType);

    if (cstrIter == phaseThermoConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName<< " type "
            << modelType << nl << nl
            << "Valid  " << typeName << " types are : " << endl
            << phaseThermoConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<cavitationModel>
    (
        cstrIter()
        (
            dict,
            alpha1,
            thermo1,
            alpha2,
            thermo2
        )
    );
}


Foam::autoPtr<Foam::cavitationModel> Foam::cavitationModel::New
(
    const dictionary& dict,
    const volScalarField& alpha1,
    const fluidThermo& thermo1,
    const volScalarField& alpha2,
    const fluidThermo& thermo2
)
{
    word modelType(dict.lookup(typeName));

    Info<< "Selecting " << typeName << ": " << modelType << endl;

    thermoConstructorTable::iterator cstrIter =
        thermoConstructorTablePtr_->find(modelType);

    if (cstrIter == thermoConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName<< " type "
            << modelType << nl << nl
            << "Valid  " << typeName << "types are : " << endl
            << thermoConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<cavitationModel>
    (
        cstrIter()
        (
            dict,
            alpha1,
            thermo1,
            alpha2,
            thermo2
        )
    );
}


// ************************************************************************* //

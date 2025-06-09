/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "surfaceTemperatureModel.H"
#include "isoThermalSurfaceTemperatureModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceTemperatureModel>
Foam::surfaceTemperatureModel::New
(
    const dictionary& dict,
    const phaseModel& phase
)
{
    word surfaceTemperatureModelType
    (
        dict.lookupOrDefault
        (
            "surfaceTemperatureModel",
            surfaceTemperatureModels::isoThermal::typeName
        )
    );

    Info<< "Selecting surfaceTemperatureModel: "
        << surfaceTemperatureModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(surfaceTemperatureModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown surfaceTemperatureModel type "
            << surfaceTemperatureModelType << endl << endl
            << "Valid surfaceTemperatureModel types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()
    (
        dict.optionalSubDict
        (
            surfaceTemperatureModelType + "SurfaceTemperatureCoeffs"
        ),
        phase
    );
}


// ************************************************************************* //

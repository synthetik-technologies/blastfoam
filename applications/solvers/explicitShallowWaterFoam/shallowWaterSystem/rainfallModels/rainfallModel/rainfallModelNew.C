/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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

#include "rainfallModel.H"
#include "constantRainfallModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::rainfallModel> Foam::rainfallModel::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    if (dict.found("R") && !dict.found("rainfallModel"))
    {
        Info<< "Selecting rainfallModel: "
            << rainfallModels::constant::typeName << endl;
        return autoPtr<rainfallModel>
        (
            new rainfallModels::constant(mesh, dict)
        );
    }

    const word model(dict.lookup("rainfallModel"));
    Info<< "Selecting rainfallModel: " << model << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(model);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown rainfallModel type "
            << model << nl << nl
            << "Valid rainfallModel for no derivatives are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return autoPtr<rainfallModel>(cstrIter()(mesh, dict));
}


// ************************************************************************* //

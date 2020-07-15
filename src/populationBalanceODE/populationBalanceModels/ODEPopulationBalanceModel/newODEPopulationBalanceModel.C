/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 Alberto Passalacqua
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

#include "ODEPopulationBalanceModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::ODEPopulationBalanceModel>
Foam::ODEPopulationBalanceModel::New
(
    const word& name,
    const dictionary& dict,
    const surfaceScalarField& phi
)
{
    word ODEPopulationBalanceModelType(dict.lookup("populationBalanceModel"));

    Info<< "Selecting ODEPopulationBalanceModel "
        << ODEPopulationBalanceModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(ODEPopulationBalanceModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown ODEPopulationBalanceModelType type "
            << ODEPopulationBalanceModelType << endl << endl
            << "Valid ODEPopulationBalanceModelType types are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return
        autoPtr<ODEPopulationBalanceModel>
        (
            cstrIter()
            (
                name,
                dict.subDict(ODEPopulationBalanceModelType + "Coeffs"),
                phi
            )
        );
}


// ************************************************************************* //

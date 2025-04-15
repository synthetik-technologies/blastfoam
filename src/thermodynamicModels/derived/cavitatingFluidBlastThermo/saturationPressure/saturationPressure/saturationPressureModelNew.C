/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
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

#include "saturationPressureModel.H"
#include "constantSaturationPressureModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::saturationPressureModel> Foam::saturationPressureModel::New
(
    const word& name,
    const dictionary& dict
)
{
    word modelTypeName;
    const dictionary* coeffDictPtr;
    if (!dict.isDict(name))
    {
        Istream& is(dict.lookup(name));
        token t(is);
        if (!t.isWord())
        {
            return autoPtr<saturationPressureModel>
            (
                new saturationPressureModels::constant
                (
                    dict,
                    t.number()
                )
            );
        }

        modelTypeName = t.wordToken();
        coeffDictPtr = &dict.optionalSubDict(modelTypeName + "Coeffs");
    }
    else
    {
        modelTypeName = dict.optionalSubDict(name).lookup<word>("type");
        coeffDictPtr = &dict.subDict(name).optionalSubDict(modelTypeName + "Coeffs");
    }

    Info<< "Selecting saturationPressureModel: " << modelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown saturationPressureModel type "
            << modelTypeName << endl << endl
            << "Valid saturationPressureModel types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(*coeffDictPtr);
}


// ************************************************************************* //

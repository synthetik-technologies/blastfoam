/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "burstModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstModel, 0);
    defineRunTimeSelectionTable(burstModel, dictionary);
}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::burstModel> Foam::burstModel::New
(
    const dictionary& dict
)
{
    word modelType(dict.lookup("burstModel"));

    Info<< "Selecting burst model: " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown burst model "
            << modelType << endl << endl
            << "Valid burst models are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burstModel::burstModel(const dictionary& dict)
:
    partialBurst_(dict.lookup<bool>("partialBurst")),
    useDelta_(dict.lookupOrDefault<bool>("useDelta", true)),
    burst_(false)
{}


Foam::burstModel::burstModel()
:
    partialBurst_(false),
    burst_(false),
    useDelta_(true)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstModel::~burstModel()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

void Foam::burstModel::writeData(Ostream& os) const
{
    writeEntry(os, "burstModel", type());
    writeEntry(os, "partialBurst", partialBurst_);
    writeEntry(os, "useDelta", useDelta_);
}


// ************************************************************************* //

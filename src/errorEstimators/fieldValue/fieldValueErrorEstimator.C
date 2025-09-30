/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2022
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

#include "fieldValueErrorEstimator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace errorEstimators
{
    defineTypeNameAndDebug(fieldValue, 0);
    addToRunTimeSelectionTable(errorEstimator, fieldValue, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorEstimators::fieldValue::fieldValue
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
:
    errorEstimator(mesh, dict, name),
    fieldName_(this->lookupFieldName(dict, typeName))
{
    this->read(dict);
}


Foam::errorEstimators::fieldValue::fieldValue
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name,
    const word& fieldName
)
:
    errorEstimator(mesh, dict, name),
    fieldName_(fieldName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimators::fieldValue::~fieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::errorEstimators::fieldValue::update(const bool scale)
{
    if (updateCurTimeIndex(!scale))
    {
        return;
    }

    const labelHashSet& eCells = this->errorCells();

    tmp<volScalarField> terrorFld(this->getFieldValue(fieldName_, eCells));
    const volScalarField& errorFld = terrorFld();

    forAllConstIter(labelHashSet, eCells, iter)
    {
        error_[iter.key()] = errorFld[iter.key()];
    }

    if (scale)
    {
        normalize(error_, eCells);
    }
}

// ************************************************************************* //

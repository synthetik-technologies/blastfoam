/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026
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

#include "linearCombinationErrorEstimator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace errorEstimators
{
    defineTypeNameAndDebug(linearCombination, 0);
    addToRunTimeSelectionTable(errorEstimator, linearCombination, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::errorEstimators::linearCombination::readCoeffs
(
    const dictionary& dict,
    wordList& fields,
    scalarList& coeffs
)
{
    fields = dict.toc();
    coeffs.setSize(fields.size());
    forAll(fields, i)
    {
        coeffs[i] = dict.lookup<scalar>(fields[i]);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::errorEstimators::linearCombination::sum
(
    const wordList& fields,
    const scalarList& coeffs
) const
{
    tmp<volScalarField> tsum
    (
        volScalarField::New
        (
            IOobject::groupName(name_, "sum"),
            mesh_,
            dimensionedScalar(dimless, 0.0)
        )
    );
    volScalarField& s = tsum.ref();

    forAll(fields, i)
    {
        if (mesh_.foundObject<volScalarField>(fields[i]))
        {
            s += coeffs[i]*mesh_.lookupObject<volScalarField>(fields[i]);
        }
    }

    return tsum;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorEstimators::linearCombination::linearCombination
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
:
    errorEstimator(mesh, dict, name),
    small_(small)
{
    this->read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimators::linearCombination::~linearCombination()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::errorEstimators::linearCombination::update(const bool scale)
{
    if (updateCurTimeIndex(!scale))
    {
        return;
    }

    const labelHashSet& eCells = this->errorCells();

    tmp<volScalarField> tnum(sum(fields_, coeffs_));
    if (denFields_.size())
    {
        tmp<volScalarField> tden(sum(denFields_, denCoeffs_));
        error_ = tnum()/max(tden(), dimensionedScalar(dimless, small_));
    }
    else
    {
        error_ = tnum();
    }

    if (scale)
    {
        normalize(error_, eCells);
    }
}


void Foam::errorEstimators::linearCombination::read(const dictionary& dict)
{
    errorEstimator::read(dict);

    readCoeffs(dict.subDict("fields"), fields_, coeffs_);
    if (dict.found("denominator"))
    {
        readCoeffs(dict.subDict("denominator"), denFields_, denCoeffs_);
    }
    small_ = dict.lookupOrDefault<scalar>("small", small);
}

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technologies
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

#include "multicomponent.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace errorEstimators
{
    defineTypeNameAndDebug(multicomponent, 0);
    addToRunTimeSelectionTable(errorEstimator, multicomponent, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorEstimators::multicomponent::multicomponent
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    errorEstimator(mesh, dict),
    names_(dict.lookup("estimators")),
    errors_(names_.size())
{
    forAll(names_, i)
    {
        errors_.set
        (
            i,
            errorEstimator::New(mesh, dict.subDict(names_[i])).ptr()
        );
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimators::multicomponent::~multicomponent()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::errorEstimators::multicomponent::read(const dictionary& dict)
{
    forAll(errors_, i)
    {
        errors_[i].read(dict.subDict(names_[i]));
    }
}


void Foam::errorEstimators::multicomponent::update()
{
    volScalarField error
    (
        IOobject
        (
            "error",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        -1.0
    );

    forAll(errors_, i)
    {
        errors_[i].update();
        error = max(error,  errors_[i].error());
    }
    error_ = error;
}

// ************************************************************************* //

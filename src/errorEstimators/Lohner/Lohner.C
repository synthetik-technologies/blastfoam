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

#include "Lohner.H"
#include "fvc.H"
#include "cubic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace errorEstimators
{
    defineTypeNameAndDebug(Lohner, 0);
    addToRunTimeSelectionTable(errorEstimator, Lohner, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorEstimators::Lohner::Lohner
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    errorEstimator(mesh, dict),
    fieldName_(dict.lookup("deltaField")),
    x_(mesh.lookupObject<volScalarField>("rho")),
    epsilon_(readScalar(dict.lookup("epsilon")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimators::Lohner::~Lohner()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::errorEstimators::Lohner::update()
{
    volScalarField& error(*this);
    surfaceScalarField xf
    (
        cubic<scalar>(x_.mesh()).interpolate(x_)
    );

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const label nInternalFaces = mesh_.nInternalFaces();
    error = 0.0;

    vector solutionD((vector(mesh_.solutionD()) + vector::one)/2.0);

    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        scalar eT =
            sqrt
            (
                mag(x_[nei] - 2.0*xf[facei] + x_[own])
               /(
                    mag(x_[nei] - xf[facei])
                  + mag(xf[facei] - x_[own])
                  + epsilon_
                   *(
                        mag(x_[nei])
                      + 2.0*mag(xf[facei])
                      + mag(x_[own])
                    )
                )
            );
        error[own] = Foam::max(error[own], eT);
        error[nei] = Foam::max(error[nei], eT);
    }
    error.correctBoundaryConditions();
}

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2020
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
    const dictionary& dict,
    const word& name
)
:
    errorEstimator(mesh, dict, name),
    fieldName_(dict.lookup("deltaField")),
    epsilon_(readScalar(dict.lookup("epsilon")))
{
    this->read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimators::Lohner::~Lohner()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::errorEstimators::Lohner::update(const bool scale)
{
    const volScalarField& x = mesh_.lookupObject<volScalarField>(fieldName_);
    surfaceScalarField xf
    (
        cubic<scalar>(mesh_).interpolate(x)
    );

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const label nInternalFaces = mesh_.nInternalFaces();
    error_ = 0.0;

    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        scalar eT =
            sqrt
            (
                mag(x[nei] - 2.0*xf[facei] + x[own])
               /(
                    mag(x[nei] - xf[facei])
                  + mag(xf[facei] - x[own])
                  + epsilon_
                   *(
                        mag(x[nei])
                      + 2.0*mag(xf[facei])
                      + mag(x[own])
                    )
                )
            );
        error_[own] = Foam::max(error_[own], eT);
        error_[nei] = Foam::max(error_[nei], eT);
    }

    forAll(error_.boundaryField(), patchi)
    {
        if (error_.boundaryField()[patchi].coupled())
        {
            const fvPatch& patch = x.boundaryField()[patchi].patch();

            const labelUList& faceCells = patch.faceCells();
            scalarField xp
            (
                x.boundaryField()[patchi].patchInternalField()
            );
            scalarField xn
            (
                x.boundaryField()[patchi].patchNeighbourField()
            );
            scalarField xbf
            (
                x.boundaryField()[patchi]
            );

            forAll(faceCells, facei)
            {
               scalar eT =
                    sqrt
                    (
                        mag(xn[facei] - 2.0*xbf[facei] + xp[facei])
                       /(
                            mag(xn[facei] - xbf[facei])
                          + mag(xbf[facei] - xp[facei])
                          + epsilon_
                           *(
                                mag(xn[facei])
                              + 2.0*mag(xbf[facei])
                              + mag(xp[facei])
                            )
                        )
                    );
                error_[faceCells[facei]] =
                    Foam::max(error_[faceCells[facei]], eT);
            }
        }
    }
    if (scale)
    {
        normalize(error_);
    }
}

// ************************************************************************* //

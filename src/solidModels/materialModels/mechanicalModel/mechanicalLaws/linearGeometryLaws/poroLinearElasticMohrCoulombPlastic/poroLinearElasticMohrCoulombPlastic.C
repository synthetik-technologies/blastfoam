/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "poroLinearElasticMohrCoulombPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "mechanicalModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(poroLinearElasticMohrCoulombPlastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, poroLinearElasticMohrCoulombPlastic, linGeomMechLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::poroLinearElasticMohrCoulombPlastic::poroLinearElasticMohrCoulombPlastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    linearElasticMohrCoulombPlastic(name, mesh, dict, nonLinGeom)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::poroLinearElasticMohrCoulombPlastic::
~poroLinearElasticMohrCoulombPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::poroLinearElasticMohrCoulombPlastic::correct
(
    volSymmTensorField& sigma
)
{
    // Call Mohr-Coulomb law to calculate the effective stress
    linearElasticMohrCoulombPlastic::correct(sigma);

    // Lookup the pore-pressure from the solver
    const volScalarField& p = mesh().lookupObject<volScalarField>("p");

    // The total stress is equal to the sum of the effective stress and
    // pore-pressure components
    sigma -= p*symmTensor(I);
}


void Foam::poroLinearElasticMohrCoulombPlastic::correct
(
    surfaceSymmTensorField& sigma
)
{
    // Call Mohr-Coulomb law to calculate the effective stress
    linearElasticMohrCoulombPlastic::correct(sigma);

    // Lookup the pressure field from the solver
    const volScalarField& p = mesh().lookupObject<volScalarField>("p");

    // Interpolate pressure to the faces
    const surfaceScalarField pf(fvc::interpolate(p));

    // The total stress is equal to the sum of the effective stress and
    // pore-pressure components
    sigma -= pf*symmTensor(I);
}


// ************************************************************************* //

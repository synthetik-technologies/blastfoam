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

#include "poroLinearElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "mechanicalModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(poroLinearElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, poroLinearElastic, linGeomMechLaw
    );
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::poroLinearElastic::makeP0f() const
{
    if (p0fPtr_)
    {
        FatalErrorIn("void Foam::poroLinearElastic::makeP0f() const")
            << "pointer already set" << abort(FatalError);
    }

    p0fPtr_ =
        new surfaceScalarField
        (
            "p0f",
            fvc::interpolate(p0_)
        );
}


const Foam::surfaceScalarField& Foam::poroLinearElastic::p0f() const
{
    if (!p0fPtr_)
    {
        makeP0f();
    }

    return *p0fPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::poroLinearElastic::poroLinearElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    linearElastic(name, mesh, dict, nonLinGeom),
    p0_
    (
        IOobject
        (
            "p0",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dict.lookupOrDefault<dimensionedScalar>
        (
            "p0",
            dimensionedScalar("zero", dimPressure, 0.0)
        )
    ),
    p0fPtr_(NULL)
{
    if (gMax(mag(p0_)()) > SMALL)
    {
        Info<< "Reading p0 initial/residual pore-pressure field" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::poroLinearElastic::~poroLinearElastic()
{
    deleteDemandDrivenData(p0fPtr_);
}


void Foam::poroLinearElastic::correct(volSymmTensorField& sigma)
{
    // Calculate effective stress
    linearElastic::correct(sigma);

    // Lookup the pressure field from the solver
    const volScalarField& p = mesh().lookupObject<volScalarField>("p");

    // Calculate the total stress as the sum of the effective stress and the
    // pore-pressure
    sigma -= (p + p0_)*symmTensor(I);
}


void Foam::poroLinearElastic::correct(surfaceSymmTensorField& sigma)
{
    // Calculate effective stress
    linearElastic::correct(sigma);

    // Lookup the pressure field from the solver
    const volScalarField& p = mesh().lookupObject<volScalarField>("p");

    // Interpolate pressure to the faces
    const surfaceScalarField pf(fvc::interpolate(p));

    // Calculate the total stress as the sum of the effective stress and the
    // pore-pressure
    sigma -= (pf + p0f())*symmTensor(I);
}


// ************************************************************************* //

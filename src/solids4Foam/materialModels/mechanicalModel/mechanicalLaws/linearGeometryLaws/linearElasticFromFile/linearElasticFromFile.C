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

#include "linearElasticFromFile.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElasticFromFile, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, linearElasticFromFile, linGeomMechLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::linearElasticFromFile::linearElasticFromFile
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    E_
    (
        IOobject
        (
            "E",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    nu_("nu", dimless, dict),
    mu_
    (
        IOobject
        (
            "mu",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        E_/(1.0 + nu_)
    ),
    lambda_
    (
        IOobject
        (
            "mu",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        E_*nu_/((1.0 + nu_)*(1.0 - 2.0*nu_))
    ),
    muf_(fvc::interpolate(mu_)),
    lambdaf_(fvc::interpolate(lambda_))
{
    if (planeStress())
    {
        lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - nu_));
        lambdaf_ = fvc::interpolate(lambda_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearElasticFromFile::~linearElasticFromFile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::linearElasticFromFile::impK() const
{
    return volScalarField::New
    (
        "impK",
        2.0*mu_ + lambda_
    );
}


Foam::tmp<Foam::scalarField>
Foam::linearElasticFromFile::impK(const label patchi) const
{
    return
        2.0*mu_.boundaryField()[patchi]
      + lambda_.boundaryField()[patchi];
}

Foam::tmp<Foam::volScalarField>
Foam::linearElasticFromFile::bulkModulus() const
{
    return volScalarField::New
    (
        "bulkModulus",
        lambda_ + 2.0/3.0*mu_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::linearElasticFromFile::elasticModulus() const
{
    return volScalarField::New
    (
        "elasticModulus",
        lambda_ + 2.0*mu_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::linearElasticFromFile::shearModulus() const
{
    return volScalarField::New
    (
        "shearModulus",
        mu_
    );
}


void Foam::linearElasticFromFile::correct(volSymmTensorField& sigma)
{
    if (incremental())
    {
        FatalErrorIn
        (
            type() + "::correct(volSymmTensorField& sigma)"
        )   << "Not implemented for incremental solid solver"
            << abort(FatalError);
    }

    // Lookup gradient of displacment from the solver
    const volTensorField& gradD =
        mesh().lookupObject<volTensorField>("grad(D)");

    // Calculate stress
    sigma = mu_*twoSymm(gradD) + lambda_*tr(gradD)*I;
}


void Foam::linearElasticFromFile::correct(surfaceSymmTensorField& sigma)
{
    if (incremental())
    {
        FatalErrorIn
        (
            type() + "::correct(surfaceSymmTensorField& sigma)"
        )   << "Not implemented for incremental solid solver"
            << abort(FatalError);
    }

    // Lookup gradient of displacment from the solver
    const surfaceTensorField& gradD =
        mesh().lookupObject<surfaceTensorField>("grad(D)f");

    // Calculate stress
    sigma = muf_*twoSymm(gradD) + lambdaf_*tr(gradD)*I;
}


// ************************************************************************* //

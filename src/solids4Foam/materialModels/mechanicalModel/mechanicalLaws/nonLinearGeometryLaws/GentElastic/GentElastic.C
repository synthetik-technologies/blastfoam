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

#include "GentElastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GentElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, GentElastic, nonLinGeomMechLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::GentElastic::GentElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_(E_/(2.0*(1.0 + nu_))),
    temperature_(dict.lookup("temperature")),
    ilim_(dict.lookup("ilim")),
    N_(dict.lookup("N")),
    Vs_(dict.lookup("Vs")),
    K_
    (
        planeStress()
      ? (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_
      : (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GentElastic::~GentElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GentElastic::impK() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            ((4.0/3.0)*mu_ + K_)
        )
    );
}


void Foam::GentElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Global Constants

    // Avocado number
    const dimensionedScalar Na
    (
        "Na",
        dimensionSet(0,0,0,0,0,0,0),
        scalar(6.022141*pow(10,9)*pow(10,9)*pow(10,5))
    );

    // Volume per solvent molecule
    const dimensionedScalar omega = Vs_/Na;

    // Boltzmann constant
    const dimensionedScalar kb
    (
        "kb",
        dimensionSet(1,2,-2,-1,0,0,0),
        scalar(1.38064852/pow(10,9)/pow(10,9)/pow(10,5))
    );

    // Crosslinking density
    const dimensionedScalar N = N_/omega;

    // Jacobian of the deformation gradient
    volScalarField J = det(F());

    // Right tensor product
    const volSymmTensorField B = symm(F() & F().T());

    // Left tensor product
    const volSymmTensorField C = symm(F().T() & F());

    // Calculate the sigma Terzaghi (Gent)
    sigma = kb*temperature_*(N/J)*((ilim_/(ilim_ - tr(B) + 3))*C - I);
}


void Foam::GentElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateFf(sigma, mu_, K_))
    {
        return;
    }

    // Global Constants

    // Avocado number
    const dimensionedScalar Na
    (
        "Na",
        dimensionSet(0,0,0,0,0,0,0),
        scalar(6.022141*pow(10,9)*pow(10,9)*pow(10,5))
    );

    // Volume per solvent molecule
    const dimensionedScalar omega = Vs_/Na;

    // Boltzmann constant
    const dimensionedScalar kb
    (
        "kb",
        dimensionSet(1,2,-2,-1,0,0,0),
        scalar(1.38064852/pow(10,9)/pow(10,9)/pow(10,5))
    );

    const dimensionedScalar N = N_/omega;

    // Jacobian of the deformation gradient
    const surfaceScalarField Jf = det(Ff());

    // Right tensor product
    const surfaceSymmTensorField Bf = symm(Ff() & Ff().T());

    // Left tensor product
    const surfaceSymmTensorField Cf = symm(Ff().T() & Ff());

    // Calculate the sigma Terzaghi (Gent)
    sigma = kb*temperature_*(N/Jf)*((ilim_/(ilim_ - tr(Bf) + 3))*Cf - I);
}


// ************************************************************************* //

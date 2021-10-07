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

#include "StVenantKirchhoffElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(StVenantKirchhoffElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, StVenantKirchhoffElastic, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::StVenantKirchhoffElastic::StVenantKirchhoffElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    lambda_("lambda", dimPressure, 0.0),
    mu_("mu", dimPressure, 0.0),
    K_("K", dimPressure, 0.0)
{
    // Read mechanical properties
    dimensionedScalar E("E", dimPressure, 0.0);
    dimensionedScalar nu("nu", dimless, 0.0);
    if
    (
        dict.found("E") && dict.found("nu")
     && !dict.found("mu") && !dict.found("K")
    )
    {
        E = dimensionedScalar(dict.lookup("E"));
        nu = dimensionedScalar(dict.lookup("nu"));

        mu_ = (E/(2.0*(1.0 + nu)));
    }
    else if
    (
        dict.found("mu") && dict.found("K")
     && !dict.found("E") && !dict.found("nu")
    )
    {
        mu_ = dimensionedScalar(dict.lookup("mu"));
        K_ = dimensionedScalar(dict.lookup("K"));

        E = 9*K_*mu_/(3*K_ + mu_);
        nu = (3*K_ - 2*mu_)/(2*(3*K_ + mu_));
    }
    else
    {
        FatalErrorIn(type())
            << "Either E and nu or mu and K should be specified"
            << abort(FatalError);
    }

    if (planeStress())
    {
        lambda_ = nu*E/((1 + nu)*(1 - nu));
    }
    else
    {
        lambda_ = nu*E/((1 + nu)*(1 - 2.0*nu));
    }

    // Set bulk modulus
    K_ = lambda_ + (2.0/3.0)*mu_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::StVenantKirchhoffElastic::~StVenantKirchhoffElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::StVenantKirchhoffElastic::impK() const
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
            2.0*mu_ + lambda_
        )
    );
}


void Foam::StVenantKirchhoffElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the right Cauchy–Green deformation tensor
    const volSymmTensorField c(symm(F().T() & F()));

    // Calculate the Green strain tensor
    const volSymmTensorField E(0.5*(c - I));

    // Calculate the 2nd Piola Kirchhoff stress
    const volSymmTensorField S(2.0*mu_*E + lambda_*tr(E)*I);

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F()));

    // Convert the 2nd Piola Kirchhoff stress to the Cauchy stress
    // sigma = (1.0/J)*symm(F() & S & F().T());
    sigma = (1.0/J)*transform(F(), S);
}


void Foam::StVenantKirchhoffElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateFf(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the right Cauchy–Green deformation tensor
    const surfaceSymmTensorField c(symm(Ff().T() & Ff()));

    // Calculate the Green strain tensor
    const surfaceSymmTensorField E(0.5*(c - I));

    // Calculate the 2nd Piola Kirchhoff stress
    const surfaceSymmTensorField S(2.0*mu_*E + lambda_*tr(E)*I);

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(Ff()));

    // Convert the 2nd Piola Kirchhoff stress to the Cauchy stress
    // sigma = (1.0/J)*symm(Ff() & S & Ff().T());
    sigma = (1.0/J)*transform(Ff(), S);
}


// ************************************************************************* //

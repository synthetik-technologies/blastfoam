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

#include "neoHookeanElastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, neoHookeanElastic, nonLinGeomMechLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::neoHookeanElastic::neoHookeanElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    E_("E", dimPressure, 0.0),
    nu_("nu", dimless, 0.0),
    mu_("mu", dimPressure, 0.0),
    K_("K", dimPressure, 0.0),
    lambda_("lambda", dimPressure, 0.0)
{
    // Read mechanical properties
    if (dict.found("E") && dict.found("nu"))
    {
        E_.readIfPresent(dict);
        nu_.readIfPresent(dict);

        mu_ = (E_/(2.0*(1.0 + nu_)));

        if (planeStress())
        {
            K_ = (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_;
        }
        else
        {
            K_ = (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_;
        }
    }
    else if (dict.found("mu") && dict.found("K"))
    {
        mu_ = dimensionedScalar("mu", dimless, dict);
        K_ = dimensionedScalar("K", dimPressure, dict);
        if (planeStress())
        {
            E_ = 4.0*K_*mu_/(K_ + mu_);
            nu_ = (K_ - mu_)/(K_ + mu_);
        }
        else
        {
            E_ = 9.0*K_*mu_/(3.0*K_ + mu_);
            nu_ = (3.0*K_ - 2.0*mu_)/(2.0*(3.0*K_ + mu_));
        }
    }
    else
    {
        FatalErrorInFunction
            << "Either E and nu or mu and K should be specified"
            << abort(FatalError);
    }

    if (planeStress())
    {
        lambda_ = K_ - mu_;
    }
    else
    {
        lambda_ = K_ - 2.0*mu_/3.0;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neoHookeanElastic::~neoHookeanElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::neoHookeanElastic::impK() const
{
    return volScalarField::New
    (
        "impK",
        mesh(),
        (4.0/3.0)*mu_ + K_ // == 2*mu + lambda
    );
}


Foam::tmp<Foam::scalarField>
Foam::neoHookeanElastic::impK(const label patchi) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            mesh().C().boundaryField()[patchi].size(),
            (4.0/3.0)*mu_.value() + K_.value()
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElastic::bulkModulus() const
{
    return volScalarField::New
    (
        "bulkModulus",
        mesh(),
        K_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElastic::elasticModulus() const
{
    return volScalarField::New
    (
        "elasticModulus",
        mesh(),
        K_ + 4.0/3.0*mu_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElastic::shearModulus() const
{
    return volScalarField::New
    (
        "shearModulus",
        mesh(),
        mu_
    );
}


void Foam::neoHookeanElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const volScalarField& J(this->J());

    // Calculate the volume preserving left Cauchy Green strain
    const volSymmTensorField bEbar(pow(J, -2.0/3.0)*symm(F() & F().T()));

    // Calculate the deviatoric stress
    const volSymmTensorField s(mu_*dev(bEbar));

    // Calculate the Cauchy stress
    sigma = (0.5*K_*(sqr(J) - 1.0)*I + s)/J;
}


void Foam::neoHookeanElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateFf(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField& J(this->Jf());

    // Calculate left Cauchy Green strain tensor with volumetric term removed
    const surfaceSymmTensorField bEbar(pow(J, -2.0/3.0)*symm(Ff() & Ff().T()));

    // Calculate deviatoric stress
    const surfaceSymmTensorField s(mu_*dev(bEbar));

    // Calculate the Cauchy stress
    sigma = (0.5*K_*(sqr(J) - 1.0)*I + s)/J;
}

Foam::tmp<Foam::volTensorField>
Foam::neoHookeanElastic::P(const volSymmTensorField& sigma) const
{
    const volScalarField& J = relative() ? this->relJ() : this->J();
    const volTensorField& F = relative() ? this->relF() : this->F();
    volTensorField H(J*T(inv(F)));
    return volTensorField::New
    (
        "P",
        mu_*pow(J, -2.0/3.0)*F
      - ((mu_/3.0)*pow(J,(-5.0/3.0))*(F && F)*H) + K_*(J-1.0)*H
    );
}


Foam::tmp<Foam::surfaceTensorField>
Foam::neoHookeanElastic::P(const surfaceSymmTensorField& sigma) const
{
    const surfaceScalarField& J = relative() ? this->relJf() : this->Jf();
    const surfaceTensorField& F = relative() ? this->relFf() : this->Ff();
    surfaceTensorField H(J*T(inv(F)));
    return surfaceTensorField::New
    (
        "P",
        mu_*pow(J, -2.0/3.0)*F
      - ((mu_/3.0)*pow(J,(-5.0/3.0))*(F && F)*H) + K_*(J-1.0)*H
    );
}


// ************************************************************************* //

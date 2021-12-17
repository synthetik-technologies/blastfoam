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

#include "linearElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, linearElastic, linGeomMechLaw
    );
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::linearElastic::calculateHydrostaticStress
(
    volScalarField& sigmaHyd,
    const volScalarField& trEpsilon
)
{
    if (solvePressureEqn_)
    {
        // Store previous iteration to allow relaxation, if needed
        sigmaHyd.storePrevIter();

        // Lookup the momentum equation inverse diagonal field
        const volScalarField* ADPtr = NULL;
        if (mesh().foundObject<volScalarField>("DEqnA"))
        {
            ADPtr = &mesh().lookupObject<volScalarField>("DEqnA");
        }
        else if (mesh().foundObject<volScalarField>("DDEqnA"))
        {
            ADPtr = &mesh().lookupObject<volScalarField>("DDEqnA");
        }
        else
        {
            FatalErrorInFunction
                << "Cannot find the DEqnA or DDEqnA field: this should be "
                << "stored in the solidModel" << abort(FatalError);
        }
        const volScalarField& AD = *ADPtr;

        // Pressure diffusivity field multiple by (4.0/3.0)*mu + K, which is
        // equivalent to (2*mu + lambda)
        // Note: we can scale this coefficient by pressureSmoothingCoeff to
        // provide greater smoothing
        const surfaceScalarField rDAf
        (
            "rDAf",
            pressureSmoothingScaleFactor_
           *fvc::interpolate
            (
                ((4.0/3.0)*mu_ + K_)/AD, "interpolate(grad(sigmaHyd))"
            )
        );

        // Solve pressure laplacian
        // Note: the the laplacian and div terms combine to produce a
        // third-order smoothing/dissipation term
        // If we only used the laplacian term then the smoothing/dissipation
        // would be second-order.
        // It would be interesting to see how this compares to the JST 2nd/4th
        // order dissipation term
        fvScalarMatrix sigmaHydEqn
        (
            fvm::Sp(1.0, sigmaHyd)
          - fvm::laplacian(rDAf, sigmaHyd, "laplacian(DA,sigmaHyd)")
          + fvc::div(rDAf*fvc::interpolate(fvc::grad(sigmaHyd)) & mesh().Sf())
         ==
            K_*trEpsilon
        );

        // Solve the pressure equation
        sigmaHydEqn.solve();

        // Relax the field
        sigmaHyd.relax();
    }
    else
    {
        // Directly calculate hydrostatic stress from displacement field
        sigmaHyd = K_*trEpsilon;
    }
}


void Foam::linearElastic::calculateHydrostaticStress
(
    surfaceScalarField& sigmaHyd,
    const surfaceScalarField& trEpsilon
)
{
    if (solvePressureEqn_)
    {
        FatalErrorInFunction
            << "'solvePressureEqn' option only implemented for volField stress "
            << "calculation" << abort(FatalError);
    }
    else
    {
        // Directly calculate hydrostatic stress from displacement field
        sigmaHyd = K_*trEpsilon;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::linearElastic::linearElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    mu_("mu", dimPressure, 0.0),
    K_("K", dimPressure, 0.0),
    E_("E", dimPressure, 0.0),
    nu_("nu", dimless, 0.0),
    lambda_("lambda", dimPressure, 0.0)
{
    // Read elastic parameters
    // The user can specify E and nu or mu and K
    if (dict.found("E") && dict.found("nu"))
    {
        // Read the Young's modulus
        E_.value() = dict.lookup<scalar>("E");

        // Read the Poisson's ratio
        nu_.value() = dict.lookup<scalar>("nu");

        // Set the shear modulus
        mu_ = E_/(2.0*(1.0 + nu_));

        // Set the bulk modulus
        if (nu_.value() < 0.5)
        {
            K_ = (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_;
        }
        else
        {
            K_.value() = GREAT;
        }
    }
    else if (dict.found("mu") && dict.found("K"))
    {
        // Read shear modulus
        mu_.value() = dict.lookup<scalar>("mu");

        // Read bulk modulus
        K_.value() = dict.lookup<scalar>("K");

        // Calculate Young's modulus
        E_ = 9.0*K_*mu_/(3.0*K_ + mu_);

        // Calculate Poisson's ratio
        nu_ = (3.0*K_ - 2.0*mu_)/(2.0*(3.0*K_ + mu_));
    }
    else
    {
        FatalErrorInFunction
            << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }

    // Set first Lame parameter
    if (nu_.value() < 0.5)
    {
        lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_));
    }
    else
    {
        lambda_.value() = GREAT;
    }

    // Check for physical Poisson's ratio
    if (nu_.value() < -1.0 || nu_.value() > 0.5)
    {
        FatalErrorInFunction
            << "Unphysical Poisson's ratio: nu should be >= -1.0 and <= 0.5"
            << abort(FatalError);
    }

    // Check for incompressibility or quasi-incompressibility
    if (nu_.value() > 0.49 && !solvePressureEqn_)
    {
        WarningInFunction
            << "Poisson's ratio is greater than 0.49: "
            << "consider setting 'solvePressureEqn' to 'yes'!" << endl;
    }

    if (solvePressureEqn_)
    {
        Info<< "    Laplacian equation will be solved for pressure" << nl
            << "    pressureSmoothingCoeff: " << pressureSmoothingScaleFactor_
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearElastic::~linearElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::linearElastic::impK() const
{
    return volScalarField::New
    (
        "impK",
        mesh(),
        nu_.value() == 0.5 ? 2.0*mu_ : 2.0*mu_ + lambda_
    );
}


Foam::tmp<Foam::scalarField>
Foam::linearElastic::impK(const label patchi) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            mesh().C().boundaryField()[patchi].size(),
            nu_.value() == 0.5
          ? 2.0*mu_.value()
          : 2.0*mu_.value() + lambda_.value()
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::linearElastic::bulkModulus() const
{
    return volScalarField::New
    (
        "bulkModulus",
        mesh(),
        K_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::linearElastic::elasticModulus() const
{
    return volScalarField::New
    (
        "elasticModulus",
        mesh(),
        lambda_ + 2.0*mu_
    );
}


Foam::tmp<Foam::volScalarField> Foam::linearElastic::shearModulus() const
{
    return volScalarField::New
    (
        "shearModulus",
        mesh(),
        mu_
    );
}


const Foam::dimensionedScalar& Foam::linearElastic::mu() const
{
    return mu_;
}


const Foam::dimensionedScalar& Foam::linearElastic::K() const
{
    return K_;
}


const Foam::dimensionedScalar& Foam::linearElastic::E() const
{
    return E_;
}


const Foam::dimensionedScalar& Foam::linearElastic::nu() const
{
    return nu_;
}


const Foam::dimensionedScalar& Foam::linearElastic::lambda() const
{
    return lambda_;
}


void Foam::linearElastic::correct(volSymmTensorField& sigma)
{
//     updateEpsilon(epsilonRef(), nu_/E_, sigma);

    // Calculate total strain
    if (incremental())
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        epsilonRef() = epsilon().oldTime() + symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        epsilonRef() = symm(gradD);
    }

    // For planeStress, correct strain in the out of plane direction
    if (planeStress())
    {
        if (mesh().solutionD()[vector::Z] > -1)
        {
            FatalErrorIn
            (
                "void Foam::linearElasticMisesPlastic::"
                "correct(volSymmTensorField& sigma)"
            )   << "For planeStress, this material law assumes the empty "
                << "direction is the Z direction!" << abort(FatalError);
        }

        epsilonRef().replace
        (
            symmTensor::ZZ,
           -(nu_/E_)
           *(sigma.component(symmTensor::XX) + sigma.component(symmTensor::YY))
        );
    }

    // Hooke's law : partitioned deviatoric and dilation form
    const volScalarField trEpsilon(tr(epsilon()));
    calculateHydrostaticStress(sigmaHydRef(), trEpsilon);
    sigma = 2.0*mu_*dev(epsilon()) + sigmaHyd()*I + sigma0();
}


void Foam::linearElastic::correct(surfaceSymmTensorField& sigma)
{
    // Calculate total strain
    updateEpsilon(epsilonfRef(), nu_/E_, sigma);

    // Hooke's law : standard form
    //sigma = 2.0*mu_*epsilonf_ + lambda_*tr(epsilonf_)*I + sigma0f();

    // Hooke's law : partitioned deviatoric and dilation form
    const surfaceScalarField trEpsilon(tr(epsilonf()));
    calculateHydrostaticStress(sigmaHydfRef(), trEpsilon);
    sigma = 2.0*mu_*dev(epsilonf()) + sigmaHydf()*I + sigma0f();
}


// ************************************************************************* //

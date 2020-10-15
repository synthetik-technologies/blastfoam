/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "kineticTheoryModel.H"
#include "mathematicalConstants.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(kineticTheoryModel, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheorySystem&
Foam::kineticTheoryModel::lookupOrConstruct
(
    const fvMesh& mesh,
    const char* name
) const
{
    if (!mesh.foundObject<kineticTheorySystem>(name))
    {
        kineticTheorySystem* ktPtr
        (
            new kineticTheorySystem(phase_.fluid())
        );

        // Transfer ownership of this object to the objectRegistry
        ktPtr->store(ktPtr);
    }

    return mesh.lookupObjectRef<kineticTheorySystem>(name);
}


Foam::kineticTheoryModel::kineticTheoryModel
(
    const phaseModel& phase,
    const dictionary& dict
)
:
    phase_(phase),
    kineticTheorySystem_(lookupOrConstruct(phase.mesh(), "kineticTheorySystem")),

    maxNut_
    (
        "maxNut",
        dimensionSet(0,2,-1,0,0),
        dict.lookupOrDefault<scalar>("maxNut", 1000)
    ),

    Theta_
    (
        IOobject
        (
            IOobject::groupName("Theta", phase.name()),
            phase.fluid().mesh().time().timeName(),
            phase.fluid().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase.fluid().mesh()
    ),

    gs0_
    (
        IOobject
        (
            IOobject::groupName("gs0", phase.name()),
            Theta_.time().timeName(),
            Theta_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Theta_.mesh(),
        dimensionedScalar(dimless, 0)
    ),
    gs0Prime_
    (
        IOobject
        (
            IOobject::groupName("gs0Prime", phase.name()),
            Theta_.time().timeName(),
            Theta_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Theta_.mesh(),
        dimensionedScalar(dimless, 0)
    ),

    lambda_
    (
        IOobject
        (
            IOobject::groupName("lambdas", phase.name()),
            Theta_.time().timeName(),
            Theta_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Theta_.mesh(),
        dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0)
    ),

    Ps_
    (
        IOobject
        (
            IOobject::groupName("p", phase.name()),
            Theta_.time().timeName(),
            Theta_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Theta_.mesh(),
        dimensionedScalar(dimPressure, 0)
    ),

    Pfric_
    (
        IOobject
        (
            IOobject::groupName("Pfric", phase.name()),
            Theta_.time().timeName(),
            Theta_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Theta_.mesh(),
        dimensionedScalar(dimPressure, 0)
    ),

    Ptot_
    (
        IOobject
        (
            IOobject::groupName("Ptot", phase.name()),
            Theta_.time().timeName(),
            Theta_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Theta_.mesh(),
        dimensionedScalar(dimPressure, 0)
    ),

    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappas", phase.name()),
            Theta_.time().timeName(),
            Theta_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Theta_.mesh(),
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0), 0)
    ),

    nut_
    (
        IOobject
        (
            IOobject::groupName("nut", phase.name()),
            Theta_.time().timeName(),
            Theta_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Theta_.mesh(),
        dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0)
    ),

    nuFric_
    (
        IOobject
        (
            IOobject::groupName("nuFric", phase.name()),
            Theta_.time().timeName(),
            Theta_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Theta_.mesh(),
        dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0)
    ),

    nuTotal_
    (
        IOobject
        (
            IOobject::groupName("nuTotal", phase.name()),
            Theta_.time().timeName(),
            Theta_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Theta_.mesh(),
        dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0)
    ),
    es_(dict.lookupType<scalar>("e"))
{
    kineticTheorySystem_.addPhase(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::~kineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField>
Foam::kineticTheoryModel::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", Theta_.group()),
                Theta_.time().timeName(),
                Theta_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - (phase_.rho()*nuTotal_)*dev(twoSymm(fvc::grad(phase_.U())))
          - ((phase_.rho()*lambda_)*fvc::div(phase_.phi()))*symmTensor::I
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::kineticTheoryModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
      - fvm::laplacian(phase_.rho()*nuTotal_, U)
      - fvc::div
        (
            (phase_.rho()*nuTotal_)*dev2(T(fvc::grad(U)))
          + ((phase_.rho()*lambda_)*fvc::div(phase_.phi()))
           *dimensioned<symmTensor>("I", dimless, symmTensor::I),
            "div(sigma." + phase_.name() + ')'
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::kineticTheoryModel::pPrime() const
{
    return
        kineticTheorySystem_.dPsdAlpha(phase_)
      + phase_*kineticTheorySystem_.frictionalPressurePrime(phase_)
      + kineticTheorySystem_.frictionalPressure();
}


void Foam::kineticTheoryModel::correct()
{
    // Local references
    volScalarField alpha(max(phase_, scalar(0)));

    tmp<volTensorField> tgradU(fvc::grad(phase_.U()));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));

    // Calculating the radial distribution function
    gs0_ = kineticTheorySystem_.gs0(phase_, phase_, true);
    gs0Prime_ = kineticTheorySystem_.gs0Prime(phase_, phase_, true);

    Ps_ = kineticTheorySystem_.Ps(phase_);
    Pfric_ = alpha*kineticTheorySystem_.frictionalPressure();
    Ptot_ = Ps_ + Pfric_;

    // 'thermal' conductivity (Table 3.3, p. 49)
    kappa_ = kineticTheorySystem_.kappa(phase_, Theta_);

    // particle viscosity (Table 3.2, p.47)
    nut_ = kineticTheorySystem_.nu(phase_, Theta_);

    volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));

    // Bulk viscosity  p. 45 (Lun et al. 1984).
//     lambda_ = (4.0/3.0)*sqr(alpha)*da*gs0_*(1.0 + es_)*ThetaSqrt/sqrtPi;
    lambda_ = kineticTheorySystem_.lambda(phase_);

    // Limit viscosity and add frictional viscosity
    nut_.min(maxNut_);
    nuFric_ = min(kineticTheorySystem_.nuFrictional(), maxNut_);
    nuTotal_ = min(maxNut_, nut_ + nuFric_);

    if (debug)
    {
        Info<< typeName << ':' << nl
            << "    max(Theta) = " << max(Theta_).value() << nl
            << "    max(nut) = " << max(nuTotal_).value() << endl;
    }
}


const Foam::phaseModel& Foam::kineticTheoryModel::phase() const
{
    return phase_;
}


const Foam::volScalarField& Foam::kineticTheoryModel::Theta() const
{
    return Theta_;
}


Foam::volScalarField& Foam::kineticTheoryModel::Theta()
{
    return Theta_;
}


const Foam::volScalarField& Foam::kineticTheoryModel::gs0() const
{
    return gs0_;
}


const Foam::volScalarField& Foam::kineticTheoryModel::gs0Prime() const
{
    return gs0Prime_;
}


// ************************************************************************* //

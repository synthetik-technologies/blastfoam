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

#include "viscousHookeanElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(viscousHookeanElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, viscousHookeanElastic, linGeomMechLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::viscousHookeanElastic::viscousHookeanElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    EInf_(dict.lookup("EInfinity")),
    E_(dict.lookup("E")),
    tau_(dict.lookup("relaxationTimes")),
    gammaInf_("gammInf", dimless, 0.0),
    gamma_(E_.size(), 0.0),
    nu_(dict.lookup("nu")),
    lambda_("lambda", dimPressure, 0.0),
    mu_("mu", dimPressure, 0.0),
    k_("k", dimPressure, 0.0),
    h_(),
    hf_(),
    s_
    (
        IOobject
        (
            "s",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    sf_
    (
        IOobject
        (
            "sf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    WilliamsLandelFerryShift_
    (
        dict.lookupOrDefault<Switch>("WilliamsLandelFerry", false)
    ),
    C1_
    (
        WilliamsLandelFerryShift_
      ? dimensionedScalar
        (
            dict.subDict("WilliamsLandelFerryCoeffs").lookup("C1")
        )
      : dimensionedScalar("C1", dimless, 0.0)
    ),
    C2_
    (
        WilliamsLandelFerryShift_
      ? dimensionedScalar
        (
            dict.subDict("WilliamsLandelFerryCoeffs").lookup("C2")
        )
      : dimensionedScalar("C2", dimTemperature, 0.0)
    ),
    Tref_
    (
        WilliamsLandelFerryShift_
      ? dimensionedScalar
        (
            dict.subDict("WilliamsLandelFerryCoeffs").lookup("Tref")
        )
      : dimensionedScalar("Tref", dimTemperature, 0.0)
    )
{
    // Check E_ and tau_ are the same length
    if (E_.size() != tau_.size())
    {
        FatalErrorIn
        (
            "Foam::viscousHookeanElastic::viscousHookeanElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "The E and relaxationTimes lists should have the same length!"
            << abort(FatalError);
    }

    // Calculate relative modulii

    const dimensionedScalar E0 =
        EInf_ + dimensionedScalar("sum(E)", dimPressure, sum(E_));

    gammaInf_ = EInf_/E0;

    forAll(gamma_, i)
    {
        gamma_[i] = E_[i]/E0.value();
    }

    // Check all the relaxation times are positive
    if (min(tau_) < SMALL)
    {
        FatalErrorIn
        (
            "Foam::viscousHookeanElastic::viscousHookeanElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "All relaxation times should be positive!"
            << abort(FatalError);
    }

    // Check all the E values are positive
    if (min(E_) < SMALL)
    {
        FatalErrorIn
        (
            "Foam::viscousHookeanElastic::viscousHookeanElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "All values of stiffness E should be positive!"
            << abort(FatalError);
    }

    // Print out the relative module

    Info<< "Relative modulii" << nl
        << "    gammaInfinity: " << gammaInf_.value() << nl;

    forAll(gamma_, i)
    {
        Info<< "    gamma[" << i << "] : " << gamma_[i] << nl;
    }

    Info<< "WilliamsLandelFerry: " << WilliamsLandelFerryShift_ << endl;

    Info<< endl;

    // Create the internal stress variables for each Maxwell model

    h_.setSize(gamma_.size());
    hf_.setSize(gamma_.size());

    forAll(h_, MaxwellModelI)
    {
        h_.set
        (
            MaxwellModelI,
            new volSymmTensorField
            (
                IOobject
                (
                    "h" + Foam::name(MaxwellModelI),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
            )
        );

        hf_.set
        (
            MaxwellModelI,
            new surfaceSymmTensorField
            (
                IOobject
                (
                    "hf" + Foam::name(MaxwellModelI),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
            )
        );

        // We need to store the old time field
        h_[MaxwellModelI].oldTime();
        hf_[MaxwellModelI].oldTime();
    }

    // Store the old time s field
    s_.oldTime();
    sf_.oldTime();

    // Set initial shear modulus
    mu_ = E0/(2.0*(1.0 + nu_));

    // Check for physical Poisson's ratio
    if (nu_.value() < -1.0 || nu_.value() > 0.5)
    {
        FatalErrorIn
        (
            "Foam::viscousHookeanElastic::viscousHookeanElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Unphysical Poisson's ratio: nu should be >= -1.0 and <= 0.5"
            << abort(FatalError);
    }

    // Check for incompressibility
    if (nu_.value() == 0.5)
    {
        Info<< "Material is incompressible: make sure to use a hybrid"
            << " approach solid model" << endl;

        // Set lambda and k to GREAT
        lambda_.value() = GREAT;
        k_.value() = GREAT;
    }
    else
    {
        if (planeStress())
        {
            lambda_ = nu_*EInf_/((1.0 + nu_)*(1.0 - nu_));
        }
        else
        {
            lambda_ = nu_*EInf_/((1.0 + nu_)*(1.0 - 2.0*nu_));
        }

        k_ = lambda_ + (2.0/3.0)*gammaInf_*mu_;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::viscousHookeanElastic::~viscousHookeanElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscousHookeanElastic::impK() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures
    scalar scaleFactor = gammaInf_.value();

    forAll(gamma_, i)
    {
        scaleFactor +=
            gamma_[i]
           *Foam::exp
            (
              - mesh().time().deltaTValue()/(2.0*tau_[i])
            );
    }

    return volScalarField::New
    (
        "impK",
        mesh(),
        nu_.value() == 0.5
      ? scaleFactor*2.0*mu_
      : scaleFactor*2.0*mu_ + lambda_
    );
}


Foam::tmp<Foam::scalarField>
Foam::viscousHookeanElastic::impK(const label patchi) const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures
    scalar scaleFactor = gammaInf_.value();

    forAll(gamma_, i)
    {
        scaleFactor +=
            gamma_[i]
           *Foam::exp
            (
              - mesh().time().deltaTValue()/(2.0*tau_[i])
            );
    }

    return tmp<scalarField>
    (
        new scalarField
        (
            mesh().C().boundaryField().size(),
            nu_.value() == 0.5
          ? scaleFactor*2.0*mu_.value()
          : scaleFactor*2.0*mu_.value() + lambda_.value()
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::viscousHookeanElastic::bulkModulus() const
{
    return volScalarField::New
    (
        "bulkModulus",
        mesh(),
        k_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::viscousHookeanElastic::elasticModulus() const
{
    return volScalarField::New
    (
        "elasticModulus",
        mesh(),
        lambda_ + 2.0*mu_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::viscousHookeanElastic::shearModulus() const
{
    return volScalarField::New
    (
        "shearModulus",
        mesh(),
        mu_
    );
}


void Foam::viscousHookeanElastic::correct(volSymmTensorField& sigma)
{
    if (incremental())
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        // Calculate deviatoric component of the strain increment
        const volSymmTensorField De(dev(symm(gradDD)));

        // Calculate deviatoric component of the initial stress, based on
        // Hooke's law
        s_ = s_.oldTime() + 2.0*mu_*De;

        // Set the volumetric component of the total stress
        sigma = tr(sigma.oldTime())*symmTensor(I) + k_*tr(gradDD)*I;
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        // Calculate deviatoric component of total strain
        const volSymmTensorField e(dev(symm(gradD)));

        // Calculate deviatoric component of the initial stress, based on
        // Hooke's law
        s_ = 2.0*mu_*e;

        // Set the volumetric component of the total stress
        sigma = k_*tr(gradD)*symmTensor(I);
    }

    // Update internal stress variables, representing stress relaxations for
    // each Maxwell model

    const scalar deltaT = mesh().time().deltaTValue();

    forAll(h_, MaxwellModelI)
    {
        volSymmTensorField& hI = h_[MaxwellModelI];
        const scalar tauI = tau_[MaxwellModelI];

        hI.storePrevIter();

        if (WilliamsLandelFerryShift_)
        {
            // Here we shift the relaxation time based on the temperature
            // field using the Williams-Landel-Ferry (WLF) approximation

            // Lookup the current temperature from the solver
            const volScalarField& T = mesh().lookupObject<volScalarField>("T");

            // Calculate the WLF shift function
            const volScalarField aT
            (
                pow(10, -C1_*(T - Tref_)/(C2_ + (T - Tref_)))
            );

            //Info<< "min(aT): " << min(aT).value() << nl
            //    << "max(aT): " << max(aT).value() << nl << endl;

            // Eqn 10.3.12 in Simo and Hughes 1998, where we shift the
            // relaxation time based on the temperature field using the
            // Williams-Landel-Ferry approximation
            hI =
                Foam::exp(-deltaT/(aT*tauI))*hI.oldTime()
              + Foam::exp(-deltaT/(2.0*aT*tauI))*(s_ - s_.oldTime());
        }
        else
        {
            // Eqn 10.3.12 in Simo and Hughes 1998
            hI =
                Foam::exp(-deltaT/tauI)*hI.oldTime()
              + Foam::exp(-deltaT/(2.0*tauI))*(s_ - s_.oldTime());
        }
    }

    // Calculate the current total stress, where the volumetric term is
    // elastic and the deviatoric term is viscoelastic
    sigma += gammaInf_*s_;

    forAll(h_, MaxwellModelI)
    {
        sigma += gamma_[MaxwellModelI]*h_[MaxwellModelI];
    }
}


void Foam::viscousHookeanElastic::correct(surfaceSymmTensorField& sigma)
{
    if (incremental())
    {
        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        // Calculate deviatoric component of the strain increment
        const surfaceSymmTensorField De(dev(symm(gradDD)));

        // Calculate deviatoric component of the initial stress, based on
        // Hooke's law
        sf_ = sf_.oldTime() + 2.0*mu_*De;

        // Set the volumetric component of the total stress
        sigma = tr(sigma.oldTime())*symmTensor(I) + k_*tr(gradDD)*I;
    }
    else
    {
        // Lookup gradient of displacement
        const surfaceTensorField& gradD =
            mesh().lookupObject<surfaceTensorField>("grad(D)f");

        // Calculate deviatoric component of total strain
        const surfaceSymmTensorField e(dev(symm(gradD)));

        // Calculate deviatoric component of the initial stress, based on
        // Hooke's law
        sf_ = 2.0*mu_*e;

        // Set the volumetric component of the total stress
        sigma = k_*tr(gradD)*symmTensor(I);
    }

    // Update internal stress variables, representing stress relaxations for
    // each Maxwell model

    const scalar deltaT = mesh().time().deltaTValue();

    forAll(hf_, MaxwellModelI)
    {
        surfaceSymmTensorField& hI = hf_[MaxwellModelI];
        const scalar tauI = tau_[MaxwellModelI];

        hI.storePrevIter();

        if (WilliamsLandelFerryShift_)
        {
            // Here we shift the relaxation time based on the temperature
            // field using the Williams-Landel-Ferry (WLF) approximation

            // Lookup the current temperature from the solver
            const surfaceScalarField& T =
                mesh().lookupObject<surfaceScalarField>("T");

            // Calculate the WLF shift function
            const surfaceScalarField aT
            (
                pow(10, -C1_*(T - Tref_)/(C2_ + (T - Tref_)))
            );

            // Eqn 10.3.12 in Simo and Hughes 1998, where we shift the
            // relaxation time based on the temperature field using the
            // Williams-Landel-Ferry approximation
            hI =
                Foam::exp(-deltaT/(aT*tauI))*hI.oldTime()
              + Foam::exp(-deltaT/(2.0*aT*tauI))*(sf_ - sf_.oldTime());
        }
        else
        {
            // Eqn 10.3.12 in Simo and Hughes 1998
            hI =
                Foam::exp(-deltaT/tauI)*hI.oldTime()
              + Foam::exp(-deltaT/(2.0*tauI))*(sf_ - sf_.oldTime());
        }
    }

    // Calculate the current total stress, where the volumetric term is
    // elastic and the deviatoric term is viscoelastic
    sigma += gammaInf_*sf_;

    forAll(hf_, MaxwellModelI)
    {
        sigma += gamma_[MaxwellModelI]*hf_[MaxwellModelI];
    }
}


Foam::scalar Foam::viscousHookeanElastic::residual()
{
    // Calculate residual based on change in internal variables

    scalar res = 0.0;

    if
    (
        mesh().time().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).foundObject<surfaceTensorField>("Ff")
    )
    {
        forAll(hf_, MaxwellModelI)
        {
            res =
                max
                (
                    res,
                    gMax
                    (
                        mag
                        (
                            hf_[MaxwellModelI].primitiveField()
                          - hf_[MaxwellModelI].prevIter().primitiveField()
                        )
                    )
                   /gMax
                    (
                        SMALL
                      + mag(hf_[MaxwellModelI].prevIter().primitiveField())
                    )
                );
        }

        return res;
    }
    else
    {
        forAll(h_, MaxwellModelI)
        {
            res =
                max
                (
                    res,
                    gMax
                    (
                        mag
                        (
                            h_[MaxwellModelI].primitiveField()
                          - h_[MaxwellModelI].prevIter().primitiveField()
                        )
                    )
                   /gMax
                    (
                        SMALL
                      + mag(h_[MaxwellModelI].prevIter().primitiveField())
                    )
                );
        }

        return res;
    }
}


// ************************************************************************* //

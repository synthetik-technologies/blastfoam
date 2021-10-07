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

#include "viscoNeoHookeanElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"
#include "fixedGradientFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(viscoNeoHookeanElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, viscoNeoHookeanElastic, nonLinGeomMechLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::viscoNeoHookeanElastic::viscoNeoHookeanElastic
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
    gammaInf_(0.0),
    gamma_(E_.size(), 0.0),
    nu_(dict.lookup("nu")),
    lambda_("lambda", dimPressure, 0.0),
    mu_(EInf_/(2.0*(1.0 + nu_))),
    mu0_(dict.lookup("mu0")),
    mu1_(dict.lookup("mu1")),
    mu2_(dict.lookup("mu2")),
    alpha_(dict.lookup("alpha")),
    beta_(dict.lookup("beta")),
    K_(dict.lookup("k")),
    h_(),
    hf_(),
    H_(),
    Hf_(),
    transformH_(),
    transformHf_(),
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
    transformNeeded_
    (
        IOobject
        (
            "transformNeeded",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimPressure, I)
    ),
    transformFbar_
    (
        IOobject
        (
            "transformFbar",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimPressure, I)
    ),
    transformNeededf_
    (
        IOobject
        (
            "transformNeededf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimPressure, I)
    ),
    transformFbarf_
    (
        IOobject
        (
            "transformFbarf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimPressure, I)
    )
{
    // Check for physical Poisson's ratio
    if (nu_.value() < -1.0 || nu_.value() > 0.5)
    {
        FatalErrorIn
        (
            "Foam::viscoNeoHookeanElastic::viscoNeoHookeanElastic\n"
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
        K_.value() = GREAT;

        if (dict.found("k"))
        {
            FatalErrorIn
            (
                "Foam::viscoNeoHookeanElastic::viscoNeoHookeanElastic\n"
                "(\n"
                "    const word& name,\n"
                "    const fvMesh& mesh,\n"
                "    const dictionary& dict\n"
                ")"
            )   << "We cannot specify k and nu independently as well as "
                << "Young's/shear modulii!" << abort(FatalError);

            // PC: what is going on here?
            K_ = dimensionedScalar(dict.lookup("k"));
        }
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

        K_ = lambda_ + (2.0/3.0)*mu_;

        if (dict.found("k"))
        {
            WarningIn
            (
                "Foam::viscoNeoHookeanElastic::viscoNeoHookeanElastic\n"
                "(\n"
                "    const word& name,\n"
                "    const fvMesh& mesh,\n"
                "    const dictionary& dict\n"
                ")"
            )   << "k is directly looked up from the dictionary; "
                << "it is not calculated from mu and lambda" << endl;

            // PC: what is going on here?
            K_ = dimensionedScalar(dict.lookup("k"));
        }
    }

    // Check E_ and tau_ are the same length
    if (E_.size() != tau_.size())
    {
        FatalErrorIn
        (
            "Foam::viscoNeoHookeanElastic::viscoNeoHookeanElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "The E and tau lists should have the same length!"
            << abort(FatalError);
    }

    // Calculate relative modulii

    const scalar E0 = EInf_.value() + sum(E_);

    gammaInf_ = EInf_.value()/E0;

    forAll(gamma_, i)
    {
        gamma_[i] = E_[i]/E0;
    }

    // Check all the relaxation times are positive
    if (min(tau_) < SMALL)
    {
        FatalErrorIn
        (
            "Foam::viscoNeoHookeanElastic::viscoNeoHookeanElastic\n"
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
            "Foam::viscoNeoHookeanElastic::viscoNeoHookeanElastic\n"
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
        << "    gammaInfinity: " << gammaInf_ << nl;

    forAll(gamma_, i)
    {
        Info<< "    gamma[" << i << "] : " << gamma_[i] << nl;
    }

    Info<< endl;

    transformH_.setSize(gamma_.size());
    transformHf_.setSize(gamma_.size());

    forAll(transformH_, MaxwellModelI)
    {
        transformH_.set
        (
            MaxwellModelI,
            new volSymmTensorField
            (
                IOobject
                (
                    "transformH" + Foam::name(MaxwellModelI),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
            )
        );

        transformHf_.set
        (
            MaxwellModelI,
            new surfaceSymmTensorField
            (
                IOobject
                (
                    "transformHf" + Foam::name(MaxwellModelI),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
            )
        );
    }


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
        h_[MaxwellModelI].storeOldTime();
        hf_[MaxwellModelI].storeOldTime();
    }


    H_.setSize(gamma_.size());
    Hf_.setSize(gamma_.size());

    forAll(H_, MaxwellModelI)
    {
        H_.set
        (
            MaxwellModelI,
            new volSymmTensorField
            (
                IOobject
                (
                    "H" + Foam::name(MaxwellModelI),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
            )
        );

        Hf_.set
        (
            MaxwellModelI,
            new surfaceSymmTensorField
            (
                IOobject
                (
                    "Hf" + Foam::name(MaxwellModelI),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
            )
        );
    }

    // Store the old time s field
    s_.storeOldTime();
    sf_.storeOldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::viscoNeoHookeanElastic::~viscoNeoHookeanElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::viscoNeoHookeanElastic::impK() const
{
    if (nu_.value() == 0.5 || mesh().foundObject<volScalarField>("p"))
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
                //K_      // Previous:  (2.0*mu_)
                // PC: this is for convergence and does/should not affect
                // the results
                2.0*mu_
            )
        );
    }
    else
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
                //K_        // Previous : (2.0*mu_ + lambda_)
                // PC: this is for convergence and does/should not affect
                // the results
                2.0*mu_ + (4.0/3.0)*K_
            )
        );
    }

}


Foam::tmp<Foam::volScalarField> Foam::viscoNeoHookeanElastic::K() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "K",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            K_
        )
    );
}


void Foam::viscoNeoHookeanElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F()));

    // Calculate the volume preserving right Cauchy Green tensor
    const volTensorField C(F().T() & F());
    // PC: this should be symmetric, we could do this:
    // const volTensorField C = symm(F().T() & F());

    // Define Fbar := J^(-1/3)*F
    const volTensorField Fbar(pow(J, -1.0/3.0)*F());
    // Define Cbar := J^(-2/3)*C
    const volTensorField Cbar(pow(J, -2.0/3.0)*C);

    // Take references to the internal fields for efficiency
    const tensorField& CI = C.primitiveField();
    const scalarField& JI = J.primitiveField();
    symmTensorField& transformNeededI = transformNeeded_.primitiveFieldRef();

    // Partial derive of WBar(CBar), internalField operation

    vector eigenValues = vector::zero;
    vector eigenValuesBar = vector::zero;

    forAll(transformNeededI, cellI)
    {
        // Calculate principal stretches: lambda^2, which are eigenvalues of C
        eigenValues = Foam::eigenValues(CI[cellI]);

        // Calculate modified principal stretches
        eigenValuesBar.x() = pow(JI[cellI], -1.0/3)*sqrt(eigenValues.x());
        eigenValuesBar.y() = pow(JI[cellI], -1.0/3)*sqrt(eigenValues.y());
        eigenValuesBar.z() = pow(JI[cellI], -1.0/3)*sqrt(eigenValues.z());

        transformNeededI[cellI].xx() =
            mu0_.value()*pow(eigenValuesBar.x(), alpha_[0] - 1)
          + mu1_.value()*pow(eigenValuesBar.x(), alpha_[1] - 1)
          + mu2_.value()*pow(eigenValuesBar.x(), alpha_[2] - 1);

        transformNeededI[cellI].yy() =
            mu0_.value()*pow(eigenValuesBar.y(), alpha_[0] - 1)
          + mu1_.value()*pow(eigenValuesBar.y(), alpha_[1] - 1)
          + mu2_.value()*pow(eigenValuesBar.y(), alpha_[2] - 1);

        transformNeededI[cellI].zz() =
            mu0_.value()*pow(eigenValuesBar.z(), alpha_[0] - 1)
          + mu1_.value()*pow(eigenValuesBar.z(), alpha_[1] - 1)
          + mu2_.value()*pow(eigenValuesBar.z(), alpha_[2] - 1);
    }

    // Partial derive of WBar(CBar), boundaryField operation
    forAll(C.boundaryField(), patchI)
    {
        // Take references to the patch fields for efficiency
        const tensorField& CP = C.boundaryField()[patchI];
        const scalarField& JP = J.boundaryField()[patchI];
        symmTensorField& transformNeededP =
            transformNeeded_.boundaryFieldRef()[patchI];

        forAll(CP, faceI)
        {
             // Calculate principal stretches: lambda^2, which are eigenvalues
             // of C
             eigenValues = Foam::eigenValues(CP[faceI]);

             // Calculate modified principal stretches
             eigenValuesBar.x() = pow(JP[faceI], -1.0/3)*sqrt(eigenValues.x());
             eigenValuesBar.y() = pow(JP[faceI], -1.0/3)*sqrt(eigenValues.y());
             eigenValuesBar.z() = pow(JP[faceI], -1.0/3)*sqrt(eigenValues.z());

             transformNeededP[faceI].xx() =
                 mu0_.value()*pow(eigenValuesBar.x(), alpha_[0] - 1)
               + mu1_.value()*pow(eigenValuesBar.x(), alpha_[1] - 1)
               + mu2_.value()*pow(eigenValuesBar.x(), alpha_[2] - 1);

             transformNeededP[faceI].yy() =
                 mu0_.value()*pow(eigenValuesBar.y(), alpha_[0] - 1)
               + mu1_.value()*pow(eigenValuesBar.y(), alpha_[1] - 1)
               + mu2_.value()*pow(eigenValuesBar.y(), alpha_[2] - 1);

             transformNeededP[faceI].zz() =
                 mu0_.value()*pow(eigenValuesBar.z(), alpha_[0] - 1)
               + mu1_.value()*pow(eigenValuesBar.z(), alpha_[1] - 1)
               + mu2_.value()*pow(eigenValuesBar.z(), alpha_[2] - 1);
        }
    }

    // 2 * FBar * (Partial Cbar of Partial WBar) * FBar.T()
    transformFbar_ = 2*transform(Fbar, transformNeeded_);

    // dev operation
    // const volSymmTensorField devFbar = dev(transformFbar_);
    sigma = dev(transformFbar_);

    // Calculate initial stress
    const volTensorField invFbar(inv(Fbar));

    s_ = transform(invFbar, sigma);

    // Update internal stress variables, representing stress relaxations for
    // each Maxwell model

    const scalar deltaT = mesh().time().deltaTValue();

    forAll(H_, MaxwellModelI)
    {
        H_[MaxwellModelI] =
            Foam::exp(-deltaT/tau_[MaxwellModelI])*h_[MaxwellModelI].oldTime()
          - Foam::exp(-deltaT/(2.0*tau_[MaxwellModelI]))*s_.oldTime();
    }

    forAll(h_, MaxwellModelI)
    {
        h_[MaxwellModelI] =
            H_[MaxwellModelI]
          + Foam::exp(-deltaT/(2.0*tau_[MaxwellModelI]))*s_;
    }

    // Calculate the current total stress, where the volumetric term is
    // elastic and the deviatoric term is viscoelastic

    scalar gRelax = gammaInf_;

    forAll(gamma_,MaxwellModelI)
    {
        gRelax +=
            gamma_[MaxwellModelI]*Foam::exp(-deltaT/(2*tau_[MaxwellModelI]));
    }

    // Calculate hydrostatic pressure, defined in (G. A. HOLZAPFEL,1996)
    const volScalarField pressure_
    (
        K_*(1.0/(beta_.value()*J))*(1 - pow(J, -beta_.value()))
    );

    sigma += J*pressure_*I + gRelax*sigma;

    forAll(H_, MaxwellModelI)
    {
        transformH_[MaxwellModelI] = transform(Fbar, H_[MaxwellModelI]);
    }

    forAll(transformH_, MaxwellModelI)
    {
        sigma += gamma_[MaxwellModelI]*dev(transformH_[MaxwellModelI]);
    }
}


void Foam::viscoNeoHookeanElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateFf(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(Ff()));

    // Calculate the volume preserving right Cauchy Green tensor
    const surfaceTensorField C(Ff().T() & Ff());
    // PC: this should be symmetric, we could do this:
    // const surfaceTensorField C = symm(F().T() & F());

    // Define Fbar := J^(-1/3)*F
    const surfaceTensorField Fbar(pow(J, -1.0/3.0)*Ff());

    // Define Cbar := J^(-2/3)*C
    const surfaceTensorField Cbar(pow(J, -2.0/3.0)*C);

    // Take references to the internal fields for efficiency
    const tensorField& CI = C.primitiveField();
    const scalarField& JI = J.primitiveField();
    symmTensorField& transformNeededI = transformNeededf_.primitiveFieldRef();

    // Partial derive of WBar(CBar), internalField operation

    vector eigenValues = vector::zero;
    vector eigenValuesBar = vector::zero;

    forAll(transformNeededI, cellI)
    {
        // Calculate principal stretches: lambda^2, which are eigenvalues of C
        eigenValues = Foam::eigenValues(CI[cellI]);

        // Calculate modified principal stretches
        eigenValuesBar.x() = pow(JI[cellI], -1.0/3)*sqrt(eigenValues.x());
        eigenValuesBar.y() = pow(JI[cellI], -1.0/3)*sqrt(eigenValues.y());
        eigenValuesBar.z() = pow(JI[cellI], -1.0/3)*sqrt(eigenValues.z());

        transformNeededI[cellI].xx() =
            mu0_.value()*pow(eigenValuesBar.x(), alpha_[0] - 1)
          + mu1_.value()*pow(eigenValuesBar.x(), alpha_[1] - 1)
          + mu2_.value()*pow(eigenValuesBar.x(), alpha_[2] - 1);

        transformNeededI[cellI].yy() =
            mu0_.value()*pow(eigenValuesBar.y(), alpha_[0] - 1)
          + mu1_.value()*pow(eigenValuesBar.y(), alpha_[1] - 1)
          + mu2_.value()*pow(eigenValuesBar.y(), alpha_[2] - 1);

        transformNeededI[cellI].zz() =
            mu0_.value()*pow(eigenValuesBar.z(), alpha_[0] - 1)
          + mu1_.value()*pow(eigenValuesBar.z(), alpha_[1] - 1)
          + mu2_.value()*pow(eigenValuesBar.z(), alpha_[2] - 1);
    }

    // Partial derive of WBar(CBar), boundaryField operation
    forAll(C.boundaryField(), patchI)
    {
        // Take references to the patch fields for efficiency
        const tensorField& CP = C.boundaryField()[patchI];
        const scalarField& JP = J.boundaryField()[patchI];
        symmTensorField& transformNeededP =
            transformNeededf_.boundaryFieldRef()[patchI];

        forAll(CP, faceI)
        {
             // Calculate principal stretches: lambda^2, which are eigenvalues
             // of C
             eigenValues = Foam::eigenValues(CP[faceI]);

             // Calculate modified principal stretches
             eigenValuesBar.x() = pow(JP[faceI], -1.0/3)*sqrt(eigenValues.x());
             eigenValuesBar.y() = pow(JP[faceI], -1.0/3)*sqrt(eigenValues.y());
             eigenValuesBar.z() = pow(JP[faceI], -1.0/3)*sqrt(eigenValues.z());

             transformNeededP[faceI].xx() =
                 mu0_.value()*pow(eigenValuesBar.x(), alpha_[0] - 1)
               + mu1_.value()*pow(eigenValuesBar.x(), alpha_[1] - 1)
               + mu2_.value()*pow(eigenValuesBar.x(), alpha_[2] - 1);

             transformNeededP[faceI].yy() =
                 mu0_.value()*pow(eigenValuesBar.y(), alpha_[0] - 1)
               + mu1_.value()*pow(eigenValuesBar.y(), alpha_[1] - 1)
               + mu2_.value()*pow(eigenValuesBar.y(), alpha_[2] - 1);

             transformNeededP[faceI].zz() =
                 mu0_.value()*pow(eigenValuesBar.z(), alpha_[0] - 1)
               + mu1_.value()*pow(eigenValuesBar.z(), alpha_[1] - 1)
               + mu2_.value()*pow(eigenValuesBar.z(), alpha_[2] - 1);
        }
    }

    // 2 * FBar * (Partial Cbar of Partial WBar) * FBar.T()
    transformFbarf_ = 2*transform(Fbar, transformNeededf_);

    // dev operation
    // const surfaceSymmTensorField devFbar = dev(transformFbar_);
    sigma = dev(transformFbarf_);

    // Calculate initial stress
    const surfaceTensorField invFbar(inv(Fbar));

    sf_ = transform(invFbar, sigma);

    // Update internal stress variables, representing stress relaxations for
    // each Maxwell model

    const scalar deltaT = mesh().time().deltaTValue();

    forAll(Hf_, MaxwellModelI)
    {
        Hf_[MaxwellModelI] =
            Foam::exp(-deltaT/tau_[MaxwellModelI])*hf_[MaxwellModelI].oldTime()
          - Foam::exp(-deltaT/(2.0*tau_[MaxwellModelI]))*sf_.oldTime();
    }

    forAll(hf_, MaxwellModelI)
    {
        hf_[MaxwellModelI] =
            Hf_[MaxwellModelI]
          + Foam::exp(-deltaT/(2.0*tau_[MaxwellModelI]))*sf_;
    }

    // Calculate the current total stress, where the volumetric term is
    // elastic and the deviatoric term is viscoelastic

    scalar gRelax = gammaInf_;

    forAll(gamma_, MaxwellModelI)
    {
        gRelax +=
            gamma_[MaxwellModelI]*Foam::exp(-deltaT/(2*tau_[MaxwellModelI]));
    }

    // Calculate hydrostatic pressure, defined in (G. A. HOLZAPFEL,1996)
    const surfaceScalarField pressure_
    (
        K_*(1.0/(beta_.value()*J))*(1 - pow(J, -beta_.value()))
    );

    sigma += J*pressure_*I + gRelax*sigma;

    forAll(Hf_, MaxwellModelI)
    {
        transformHf_[MaxwellModelI] = transform(Fbar, Hf_[MaxwellModelI]);
    }

    forAll(transformHf_, MaxwellModelI)
    {
        sigma += gamma_[MaxwellModelI]*dev(transformHf_[MaxwellModelI]);
    }
}


// ************************************************************************* //

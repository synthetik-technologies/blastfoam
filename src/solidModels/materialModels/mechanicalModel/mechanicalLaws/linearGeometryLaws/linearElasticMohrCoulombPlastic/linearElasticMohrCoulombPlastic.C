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

#include "linearElasticMohrCoulombPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElasticMohrCoulombPlastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, linearElasticMohrCoulombPlastic, linGeomMechLaw
    );
}


// * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * * //

void Foam::linearElasticMohrCoulombPlastic::calculateStress
(
    symmTensor& sigma,
    scalar& activeYield
) const
{
    // Principal stresses
    vector sigma_prin = vector::zero;

    // Principal directions
    tensor ev = tensor::zero;

    activeYield = 0;

    for (int i = 0; i < 6; i++)
    {
        if (sigma[i] > SMALL)
        {
            sigma[i] += (i + 1)*(1e-6)*sigma[i];
        }
    }

    calculateEigens(sigma_prin, ev, sigma);

    // Re-order the principal stresses
    scalar sigma1 = sigma_prin[0];
    scalar sigma2 = sigma_prin[1];
    scalar sigma3 = sigma_prin[2];

    label sigma1_po = 0;
    label sigma2_po = 1;
    label sigma3_po = 2;

    for (int i = 1; i < 3; i++)
    {
        if (sigma_prin[i] > sigma1)
        {
            sigma1 = sigma_prin[i];
            sigma1_po = i;
        }
    }

    for (int i = 1; i >= 0; i--)
    {
        if (sigma_prin[i] < sigma3)
        {
            sigma3 = sigma_prin[i];
            sigma3_po = i;
        }
    }

    for (int i = 0; i < 3; i++)
    {
        if ((i != sigma1_po) && (i != sigma3_po))
        {
            sigma2_po = i;
            sigma2 = sigma_prin[i];
        }
    }

    vector sigmaB = vector(sigma1, sigma2, sigma3);

    // Evaluate the yielding function f
    const scalar f = k_*sigmaB[0] - sigmaB[2] - 2*c_.value()*Foam::sqrt(k_);

    // Check if yielding
    if (f > SMALL)
    {
        // Determine the return type

        // Calculate the boundary planes
        const scalar t1 =
            (r_lg1_ & (invC_ & (sigma_prin - sigma_a_)))
            /(r_lg1_ & (invC_ & r_lf1_));

        const scalar t2 =
            (r_lg2_ & (invC_ & (sigma_prin - sigma_a_)))
            /(r_lg2_ & (invC_ & r_lf2_));

        const scalar p_12 = (rp_^r_lf1_) & (sigma_prin - sigma_a_);

        const scalar p_13 = (rp_^r_lf2_) & (sigma_prin - sigma_a_);

        // Compute the stress field based on the correct return type

        if ((p_12 >= 0) && (p_13 <= 0))
        {
            sigmaB -= f*rp_;
        }
        else if ((p_12 < 0) && (p_13 < 0))
        {
            sigmaB = t1*r_lf1_ + sigma_a_;
        }
        else if ((p_12 > 0) && (p_13 > 0))
        {
            sigmaB = t2*r_lf2_ + sigma_a_;
        }
        else if ((t1 > 0) && (t2 > 0))
        {
            sigmaB = sigma_a_;
        }

        sigma_prin[sigma1_po] = sigmaB[0];
        sigma_prin[sigma2_po] = sigmaB[1];
        sigma_prin[sigma3_po] = sigmaB[2];

        // Form the diagonal stress
        const diagTensor sigma_diag =
            diagTensor(sigma_prin[0], sigma_prin[1], sigma_prin[2]);

        // Transform the returned stress back into general space
        sigma = symm(ev.T() & sigma_diag & ev);

        activeYield = 1.0;
    }
}


void Foam::linearElasticMohrCoulombPlastic::calculateEigens
(
    vector& sigma_prin,
    tensor& ev,
    const symmTensor sigma
) const
{
    // Check for a zero stress tensor
    if (mag(sigma) < SMALL)
    {
        ev = I;

        return;
    }

    scalar i = 0;
    scalar ii = 0;
    scalar iii = 0;

    if
    (
        (
            mag(sigma.xy()) + mag(sigma.xz()) + mag(sigma.xy())
          + mag(sigma.yz()) + mag(sigma.xz()) + mag(sigma.yz())
        ) < 1e-6*(mag(sigma.xx()) + mag(sigma.yy()) + mag(sigma.zz()))
        )
    {
        // Diagonal matrix
        i = sigma.xx();
        ii = sigma.yy();
        iii = sigma.zz();
    }
    else
    {
        const scalar a = -sigma.xx() - sigma.yy() - sigma.zz();

        const scalar b =
            sigma.xx()*sigma.yy()
          + sigma.xx()*sigma.zz()
          + sigma.yy()*sigma.zz()
          - sigma.xy()*sigma.xy()
          - sigma.xz()*sigma.xz()
          - sigma.yz()*sigma.yz();

        const scalar c =
          - sigma.xx()*sigma.yy()*sigma.zz()
          - sigma.xy()*sigma.yz()*sigma.xz()
          - sigma.xz()*sigma.xy()*sigma.yz()
          + sigma.xz()*sigma.yy()*sigma.xz()
          + sigma.xy()*sigma.xy()*sigma.zz()
          + sigma.xx()*sigma.yz()*sigma.yz();

        // If there is a zero root
        if (mag(c) < SMALL)
        {
            const scalar disc = Foam::max(sqr(a) - 4*b, 0.0);

            WarningIn("poroMohrCoulob::calculateEigens(...)")
                << "Stress tensor has a zero root!" << endl;

            const scalar q = -0.5*Foam::sqrt(max(scalar(0), disc));

            i = 0;
            ii = -0.5*a + q;
            iii = -0.5*a - q;
        }
        else
        {
            const scalar Q = (a*a - 3.0*b)/9.0;
            const scalar R = (2.0*a*a*a - 9.0*a*b + 27.0*c)/54.0;

            const scalar R2 = sqr(R);
            const scalar Q3 = pow3(Q);

            // Three different real roots
            if (R2 < Q3)
            {
                const scalar sqrtQ = Foam::sqrt(Q);
                const scalar theta = Foam::acos(R/(Q*sqrtQ));

                const scalar m2SqrtQ = -2*sqrtQ;
                const scalar aBy3 = a/3;

                i = m2SqrtQ*Foam::cos(theta/3) - aBy3;
                ii =
                    m2SqrtQ*Foam::cos((theta + constant::mathematical::twoPi)/3.0)
                  - aBy3;
                iii =
                    m2SqrtQ*Foam::cos((theta - constant::mathematical::twoPi)/3.0)
                  - aBy3;
            }
            else
            {
                const scalar A = Foam::cbrt(R + Foam::sqrt(R2 - Q3));

                // Three equal real roots
                if (A < SMALL)
                {
                    const scalar root = -a/3;
                    i = root;
                    ii = root;
                    iii = root;
                }
                else
                {
                    // Complex roots
                    WarningIn("poroMohrCoulob::calculateEigens(...)")
                        << "Complex eigenvalues detected for symmTensor: "
                        << sigma << nl << "Setting roots to zero!" << endl;

                    i = 0;
                    ii = 0;
                    iii = 0;
                }
            }
        }
    }

    // Sort the eigenvalues into ascending order

    if (mag(i) > mag(ii))
    {
        Swap(i, ii);
    }

    if (mag(ii) > mag(iii))
    {
        Swap(ii, iii);
    }

    if (mag(i) > mag(ii))
    {
        Swap(i, ii);
    }

    sigma_prin[0] = i;
    sigma_prin[1] = ii;
    sigma_prin[2] = iii;

    for (int j = 0; j < 3; j++)
    {
        if (mag(sigma_prin[j]) < SMALL)
        {
            if (j == 0)
            {
                ev[j*3 + 0] = 1;
                ev[j*3 + 1] = 0;
                ev[j*3 + 2] = 0;
            }
            else if (j == 1)
            {
                ev[j*3 + 0] = 0;
                ev[j*3 + 1] = 1;
                ev[j*3 + 2] = 0;
            }
            else
            {
                ev[j*3 + 0] = 0;
                ev[j*3 + 1] = 0;
                ev[j*3 + 2] = 1;
            }
        }
        else
        {
            const symmTensor A = symmTensor(sigma - sigma_prin[j]*I);

            // Calculate the sub-determinants of the 3 components
            const scalar sd0 = A.yy()*A.zz() - A.yz()*A.yz();
            const scalar sd1 = A.xx()*A.zz() - A.xz()*A.xz();
            const scalar sd2 = A.xx()*A.yy() - A.xy()*A.xy();

            const scalar magSd0 = mag(sd0);
            const scalar magSd1 = mag(sd1);
            const scalar magSd2 = mag(sd2);

            // Evaluate the eigenvector using the largest sub-determinant
            if ((magSd0 > magSd1) && (magSd0 > magSd2) && (magSd0 > SMALL))
            {
                vector newEv =
                    vector
                    (
                        1,
                        (A.yz()*A.xz() - A.zz()*A.xy())/sd0,
                        (A.yz()*A.xy() - A.yy()*A.xz())/sd0
                    );
                newEv /= mag(newEv);

                ev[j*3] = newEv[0];
                ev[j*3+1] = newEv[1];
                ev[j*3+2] = newEv[2];
            }
            else if ((magSd1 > magSd2) && (magSd1 > SMALL))
            {
                vector newEv =
                    vector
                    (
                        (A.xz()*A.yz() - A.zz()*A.xy())/sd1,
                        1,
                        (A.xz()*A.xy() - A.xx()*A.yz())/sd1
                    );
                newEv /= mag(newEv);

                ev[j*3] = newEv[0];
                ev[j*3 + 1] = newEv[1];
                ev[j*3 + 2] = newEv[2];
            }
            else if (magSd2 > SMALL)
            {
                vector newEv =
                    vector
                    (
                        (A.xy()*A.yz() - A.yy()*A.xz())/sd2,
                        (A.xy()*A.xz() - A.xx()*A.yz())/sd2,
                        1
                    );
                newEv /= mag(newEv);

                ev[j*3] = newEv[0];
                ev[j*3 + 1] = newEv[1];
                ev[j*3 + 2] = newEv[2];
            }
            else
            {
                WarningIn("poroMohrCoulob::calculateEigens(...)")
                    << "Strange things detected for stress tensor: "
                    << sigma << nl
                    << "Setting eigenvectors to (0, 0, 0)" << endl;

                ev[j*3] = 0;
                ev[j*3 + 1] = 0;
                ev[j*3 + 2] = 0;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::linearElasticMohrCoulombPlastic::linearElasticMohrCoulombPlastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    E_("E", dimPressure, dict),
    nu_("nu", dimless, dict),
    lambda_
    (
        planeStress()
      ? nu_*E_/((1.0 + nu_)*(1.0 - nu_))
      : nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))
    ),
    mu_(E_/(2.0*(1.0 + nu_))),
    K_(lambda_ + (2.0/3.0)*mu_),
    varPhi_("frictionAngle", dimless, dict),
    c_(dict.lookup("cohesion")),
    varPsi_("dilationAngle", dimless, dict),
    k_
    (
        (
            (1 + sin(varPhi_/180.0*constant::mathematical::pi))
           /(1 - sin(varPhi_/180.0*constant::mathematical::pi))
        ).value()
    ),
    m_
    (
        (
            (1 + sin(varPsi_/180.0*constant::mathematical::pi))
           /(1 - sin(varPsi_/180.0*constant::mathematical::pi))
        ).value()
    ),
    a_(k_, 0, -1),
    b_(m_, 0, -1),
    C_
    (
        (K_ + (4.0/3.0)*mu_).value(),
        (K_ - (2.0/3.0)*mu_).value(),
        (K_ - (2.0/3.0)*mu_).value(),
        (K_ + (4.0/3.0)*mu_).value(),
        (K_ - (2.0/3.0)*mu_).value(),
        (K_ + (4.0/3.0)*mu_).value()
    ),
    invC_(inv(C_)),
    rp_((C_ & b_)/(b_ & (C_ & a_))),
    r_lf1_(1, 1, k_),
    r_lf2_(1, k_, k_),
    r_lg1_(1, 1, m_),
    r_lg2_(1, m_, m_),
    sigma_a_((2.0*c_*Foam::sqrt(k_)/(k_ - 1.0)).value()*vector::one),
    sigmaEff_
    (
        IOobject
        (
            "sigmaEff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    sigmaEfff_
    (
        IOobject
        (
            "sigmaEfff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    DEpsilonP_
    (
        IOobject
        (
            "DEpsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DEpsilonPf_
    (
        IOobject
        (
            "DEpsilonPf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    epsilonP_
    (
        IOobject
        (
            "epsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    epsilonPEq_
    (
        IOobject
        (
            "epsilonPEq",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    activeYield_
    (
        IOobject
        (
            "activeYield",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearElasticMohrCoulombPlastic::~linearElasticMohrCoulombPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::linearElasticMohrCoulombPlastic::impK() const
{
    return volScalarField::New
    (
        "impK",
        mesh(),
        nu_.value() == 0.5
      ? 2.0*mu_
      : 2.0*mu_ + lambda_
    );
}


Foam::tmp<Foam::scalarField>
Foam::linearElasticMohrCoulombPlastic::impK(const label patchi) const
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


Foam::tmp<Foam::volScalarField>
Foam::linearElasticMohrCoulombPlastic::bulkModulus() const
{
    return volScalarField::New
    (
        "bulkModulus",
        mesh(),
        K_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::linearElasticMohrCoulombPlastic::elasticModulus() const
{
    return volScalarField::New
    (
        "elasticModulus",
        mesh(),
        lambda_ + 2.0*mu_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::linearElasticMohrCoulombPlastic::shearModulus() const
{
    return volScalarField::New
    (
        "shearModulus",
        mesh(),
        mu_
    );
}


const Foam::dimensionedScalar& Foam::linearElasticMohrCoulombPlastic::mu() const
{
    return mu_;
}


const Foam::dimensionedScalar&
Foam::linearElasticMohrCoulombPlastic::lambda() const
{
    return lambda_;
}


void Foam::linearElasticMohrCoulombPlastic::correct(volSymmTensorField& sigma)
{
    // Reset sigma to sigmaEff, in case a derived poro-elastic based law is
    // being used
    sigma = 1.0*sigmaEff_;

    // Store the previous iteration of sigma to allow under-relaxation and also
    // to calculate the material residual
    sigma.storePrevIter();
    sigmaEff_.storePrevIter();

    // Calculate the increment of total strain
    volSymmTensorField DEpsilon
    (
        IOobject
        (
            "DEpsilon",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    if (incremental())
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        DEpsilon = symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        // Calculate gradient of displacement increment
        const volTensorField gradDD(gradD - gradD.oldTime());

        DEpsilon = symm(gradDD);
    }

    // Calculate trial elastic stress assuming Hooke's elastic law
    const volSymmTensorField DSigmaTrial
    (
        2.0*mu_*DEpsilon + lambda_*I*tr(DEpsilon)
    );

    // Set sigma to sigma effective trial
    sigma = sigmaEff_.oldTime() + DSigmaTrial;

    // Take a reference to internal fields for efficiency
    symmTensorField& sigmaI = sigma.primitiveFieldRef();
    scalarField& activeYieldI = activeYield_.primitiveFieldRef();

    // Correct sigma internal field
    forAll(sigmaI, cellI)
    {
        calculateStress(sigmaI[cellI], activeYieldI[cellI]);
    }

    // Correct sigma on the boundary patches
    forAll (sigma.boundaryField(), patchI)
    {
        // Take references to the boundary patches for efficiency
        symmTensorField& sigmaP = sigma.boundaryFieldRef()[patchI];
        scalarField& activeYieldP = activeYield_.boundaryFieldRef()[patchI];

        if (!sigma.boundaryField()[patchI].coupled())
        {
            forAll(sigmaP, faceI)
            {
                calculateStress(sigmaP[faceI], activeYieldP[faceI]);
            }
        }
    }

    // Correct the coupled boundary conditions
    sigma.correctBoundaryConditions();
    activeYield_.correctBoundaryConditions();

    // Under-relax the stress
    sigma.relax();

    // Update the effective stress tensor
    // Note: if a poro-elasto-plastic law is used then the pore-pressure term
    // will be added after this
    sigmaEff_ = 1.0*sigma;
}


void Foam::linearElasticMohrCoulombPlastic::correct
(
    surfaceSymmTensorField& sigma
)
{
    // Reset sigma to sigmaEff, in case a derived poro-elastic based law is
    // being used
    sigma = 1.0*sigmaEfff_;

    // Store the previous iteration of sigma to allow under-relaxation and also
    // to calculate the material residual
    sigma.storePrevIter();
    sigmaEfff_.storePrevIter();

    // Calculate the increment of total strain
    surfaceSymmTensorField DEpsilon
    (
        IOobject
        (
            "DEpsilon",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    if (incremental())
    {
        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        DEpsilon = symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const surfaceTensorField& gradD =
            mesh().lookupObject<surfaceTensorField>("grad(D)f");

        // Calculate gradient of displacement increment
        const surfaceTensorField gradDD(gradD - gradD.oldTime());

        DEpsilon = symm(gradDD);
    }

    // Calculate trial elastic stress assuming Hooke's elastic law
    const surfaceSymmTensorField DSigmaTrial
    (
        2.0*mu_*DEpsilon + lambda_*I*tr(DEpsilon)
    );

    // Set sigma to sigma effective trial
    sigma = sigmaEfff_.oldTime() + DSigmaTrial;

    // Take a reference to internal fields for efficiency
    symmTensorField& sigmaI = sigma.primitiveFieldRef();
    scalarField& activeYieldI = activeYield_.primitiveFieldRef();

    const labelList& faceOwner = mesh().faceOwner();
    const labelList& faceNeighbour = mesh().faceNeighbour();

    // Correct sigma internal field
    forAll(sigmaI, faceI)
    {
        calculateStress(sigmaI[faceI], activeYieldI[faceOwner[faceI]]);

        // Update the neighbour activeYield as it is a vol field
        activeYieldI[faceOwner[faceI]] = activeYieldI[faceNeighbour[faceI]];
    }

    // Correct sigma on the boundary patches
    forAll (sigma.boundaryField(), patchI)
    {
        // Take references to the boundary patches for efficiency
        symmTensorField& sigmaP = sigma.boundaryFieldRef()[patchI];
        scalarField& activeYieldP = activeYield_.boundaryFieldRef()[patchI];

        // We calculate on surfaceField coupled boundaries
        //if (!sigma.boundaryField()[patchI].coupled())
        {
            forAll(sigmaP, faceI)
            {
                calculateStress(sigmaP[faceI], activeYieldP[faceI]);
            }
        }
    }

    // Correct the coupled boundary conditions
    activeYield_.correctBoundaryConditions();

    // Under-relax the stress
    sigma.relax();

    // Update the effective stress tensor
    // Note: if a poro-elasto-plastic law is used then the pore-pressure term
    // will be added after this
    sigmaEfff_ = 1.0*sigma;
}


Foam::scalar Foam::linearElasticMohrCoulombPlastic::residual()
{
    // Calculate residual based on change in plastic strain increment
    if
    (
        mesh().foundObject<surfaceVectorField>("grad(D)f")
     || mesh().foundObject<surfaceVectorField>("grad(DD)f")
    )
    {
        return
            gMax
            (
                mag
                (
                    sigmaEfff_.primitiveField()
                  - sigmaEfff_.prevIter().primitiveField()
                )
            )/gMax(SMALL + mag(sigmaEfff_.primitiveField()));
    }
    else
    {
        return
            gMax
            (
                mag
                (
                    sigmaEff_.primitiveField()
                  - sigmaEff_.prevIter().primitiveField()
                )
            )/gMax(SMALL + mag(sigmaEff_.primitiveField()));
    }
}


void Foam::linearElasticMohrCoulombPlastic::updateTotalFields()
{
    Info<< type() << ": updating total fields" << endl;

    // Calculate increment of plastic strain
    // This is complicated because DEpsilonP also has a volumetric term
    symmTensorField& DEpsilonPI = DEpsilonP_.primitiveFieldRef();

    // Calculate the increment of total strain
    volSymmTensorField DEpsilon
    (
        IOobject
        (
            "DEpsilon",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    if (incremental())
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        DEpsilon = symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        // Calculate gradient of displacement increment
        const volTensorField gradDD(gradD - gradD.oldTime());

        DEpsilon = symm(gradDD);
    }

    const symmTensorField& DEpsilonI = DEpsilon.primitiveField();
    const symmTensorField& sigmaEffI = sigmaEff_.primitiveField();

    const scalar mu = mu_.value();
    const scalar lambda = lambda_.value();

    forAll(DEpsilonPI, cellI)
    {
        // Off-diagonal strains
        DEpsilonPI[cellI].xy() =
            DEpsilonI[cellI].xy() - (sigmaEffI[cellI].xy()/(2.0*mu));
        DEpsilonPI[cellI].xz() =
            DEpsilonI[cellI].xz() - (sigmaEffI[cellI].xz()/(2.0*mu));
        DEpsilonPI[cellI].yz() =
            DEpsilonI[cellI].yz() - (sigmaEffI[cellI].yz()/(2.0*mu));

        // Solve a linear system (Ax = b) to calculate the strains on the
        // diagonal strains
        const tensor A =
            tensor
            (
              - (2.0*mu + lambda), - lambda, - lambda,
              - lambda, - (2.0*mu + lambda), - lambda,
              - lambda,  - lambda, - (2.0*mu + lambda)
            );

        const vector b =
            vector
            (
                sigmaEffI[cellI].xx()
             - (2.0*mu + lambda)*DEpsilonI[cellI].xx()
             - lambda*DEpsilonI[cellI].yy()
             - lambda*DEpsilonI[cellI].zz(),

                sigmaEffI[cellI].yy()
             - lambda*DEpsilonI[cellI].xx()
             - (2.0*mu + lambda)*DEpsilonI[cellI].yy()
             - lambda*DEpsilonI[cellI].zz(),

                sigmaEffI[cellI].zz()
             - lambda*DEpsilonI[cellI].xx()
             - lambda*DEpsilonI[cellI].yy()
             - (2.0*mu + lambda)*DEpsilonI[cellI].zz()
            );

        const vector x = (inv(A) & b);

        DEpsilonPI[cellI].xx() = x.x();
        DEpsilonPI[cellI].yy() = x.y();
        DEpsilonPI[cellI].zz() = x.z();
    }

    //     DEpsilon
    //   - (1.0/3.0)*I*tr(sigma - sigmaEff_.oldTime())/(2.0*mu_ + 3.0*lambda_)
    //   - dev(sigma - sigmaEff_.oldTime())/(2.0*mu_);

    // Calculate the increment of equivalent plastic strain
    const volScalarField DEpsilonPEq(sqrt((2.0/3.0)*magSqr(dev(DEpsilonP_))));

    // Update the total plastic strain
    epsilonP_ = epsilonP_.oldTime() + DEpsilonP_;

    // Update the total equivalent plastic strain
    epsilonPEq_ = epsilonPEq_.oldTime() + DEpsilonPEq;

    Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq) << endl;
    Info<< "    Max epsilonPEq is " << gMax(epsilonPEq_) << endl;

    // Count cells actively yielding
    int numCellsYielding = 0;

    forAll(activeYield_.primitiveField(), cellI)
    {
        if (activeYield_.primitiveField()[cellI] > SMALL)
        {
            numCellsYielding++;
        }
    }

    reduce(numCellsYielding, sumOp<int>());

    Info<< "    " << numCellsYielding << " cells are actively yielding"
        << nl << endl;
}


// ************************************************************************* //

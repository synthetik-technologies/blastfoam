/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "plasticityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(plasticityModel, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

plasticityModel::plasticityModel
(
    const volTensorField& F,
    const dictionary& dict
)
:
    b_ ("b", F & F.T()),

    P_
    (
        IOobject
        (
            "P",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        F.mesh(),
        dimensionedTensor("P", dimensionSet(1,-1,-2,0,0,0,0), tensor::zero)
    ),

    tau_ (F),

    vMises_
    (
        IOobject
        (
            "vMises",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        F.mesh(),
        dimensionedScalar("vMises", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    ),

    strain_p_
    (
        IOobject
        (
            "strain_p",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        F.mesh(),
        dimensionedScalar("strain_p", dimless, 0.0)
    ),

    CpInv_ (F),

    p_
    (
        IOobject
        (
            "p",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        F.mesh(),
        dimensionedScalar("p", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    ),

    energyAlgorithm_
    (
        IOobject
        (
            "energyAlgorithm",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        F.mesh(),
        dimensionedScalar
        (
            "energyAlgorithm",
            dimensionSet(1,-1,-2,0,0,0,0),
            0.0
        )
    ),

    model_(dict.lookup("plasticityModel")),

    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_(E_/(2.0*(1.0 + nu_))),
    lambda_(nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))),
    kappa_(lambda_ + (2.0/3.0)*mu_),

    Hm_(dict.lookup("Hm")),
    Ys0_(dict.lookup("Ys")),

    Ys_
    (
        IOobject
        (
            "Ys",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        F.mesh(),
        Ys0_
    ),

    pWaveSpeed_(sqrt((lambda_ + 2.0*mu_)/rho_)),
    sWaveSpeed_(sqrt(mu_/rho_))
{
    correct();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

plasticityModel::~plasticityModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void plasticityModel::correct()
{
    const fvMesh& mesh_ = P_.mesh();
    const volTensorField& F = mesh_.lookupObject<volTensorField>("F");
    const volScalarField& J = mesh_.lookupObject<volScalarField>("J");
    const volVectorField& rhoU = mesh_.lookupObject<volVectorField>("rhoU");

    energyAlgorithm_ =
        0.5*mu_*(pow(J,(-2.0/3.0))*(F && F) - 3.0)
      + (0.5*kappa_*(J - 1.0)*(J - 1.0)) + (0.5*(rhoU & rhoU)/rho_);

    operations op(mesh_);

    if (model_ == "vonMisesPlasticity")
    {
        const volTensorField Finv(inv(F));

        p_ = kappa_*(log(J)/J);
        b_ = F & CpInv_ & F.T();

        forAll(mesh_.cells(), cellID)
        {
            op.eigenStructure(b_[cellID]);
            const vector& eVal = op.eigenValue();
            const tensor& eVec = op.eigenVector();

            vector eStretchT(sqrt(eVal.x()), sqrt(eVal.y()), sqrt(eVal.z()));

            // Principle Trial Deviatoric Kirchoff Stress Vector
            vector tauDevT = vector::zero;

            for (int i=0; i<3; i++)
            {
                tauDevT[i] =
                    2.0*mu_.value()*log(eStretchT[i])
                  - (2.0/3.0)*mu_.value()*log(J[cellID]);
            }

            // Yield Criterion
            vector directionV = vector::zero;
            scalar plasticM = 0.0;

            scalar f =
                sqrt((3.0/2.0)*(tauDevT&tauDevT))
              - (Ys0_.value() + Hm_.value()*strain_p_[cellID]);

            vector tauDev = tauDevT;
            vector eStretch = vector::zero;

            // Condition for Plastic Deformation
            if (f > 0.0)
            {
                directionV = tauDevT/(sqrt(2.0/3.0)*sqrt(tauDevT & tauDevT));
                plasticM = f/(3.0*mu_.value() + Hm_.value());

                // Elastic Stretch Vector
                for (int i=0; i<3; i++)
                {
                    eStretch[i] =
                        exp(log(eStretchT[i]) - plasticM*directionV[i]);
                }

                // Principle Deviatoric Kirchoff Stress Tensor
                for (int i=0; i<3; i++)
                {
                    tauDev[i] =
                        (1.0
                      - ((2.0*mu_.value()*plasticM)
                       /(sqrt(2.0/3.0)*sqrt(tauDevT & tauDevT))))
                       *tauDevT[i];
                }

                // Update Left Cauchy Green Strain Tensor
                b_[cellID] = tensor::zero;

                for (int i=0; i<3; i++)
                {
                    b_[cellID] +=
                        (eStretch[i]*eStretch[i])
                       *( vector(eVec[3*i], eVec[3*i+1], eVec[3*i+2])
                       *vector(eVec[3*i], eVec[3*i+1], eVec[3*i+2]) );
                }

                // Update Plastic Strain
                strain_p_[cellID] += plasticM;
            }

            // Kirchoff Stress Tensor
            tau_[cellID] = tensor::zero;
            for (int i=0; i<3; i++)
            {
                tau_[cellID] +=
                    (tauDev[i] + (J[cellID]*p_[cellID]))
                   *(vector(eVec[3*i], eVec[3*i+1], eVec[3*i+2])
                   *vector(eVec[3*i], eVec[3*i+1], eVec[3*i+2]));
            }

            // Update PK1
            P_[cellID] = tau_[cellID] & Finv[cellID].T();

            // Update von-Mises stresses
            vMises_[cellID] = sqrt(1.5*(tauDev && tauDev));
        }

        CpInv_ = Finv & b_ & Finv.T();
    }

    else
    {
        FatalErrorIn
        (
            "plasticityModel.C"
        )   << "Valid type entry is 'vonMisesPlasticity' for constitutiveModel"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void plasticityModel::printMaterialProperties()
{
    Info<< "\nPrinting material properties ..." << nl
        << "Constitutive model = " << model_ << nl
        << "Density = " << rho_.value() << " " << rho_.dimensions() << nl
        << "Young's modulus = " << E_.value() << " " << E_.dimensions() << nl
        << "Poisson's ratio = " << nu_.value() << " " << nu_.dimensions() << nl
        << "Lame's first parameter lambda = " << lambda_.value() << " "
        << lambda_.dimensions() << nl
        << "Lame's second parameter mu = " << mu_.value() << " "
        << mu_.dimensions() << nl
        << "Bulk modulus kappa = " << kappa_.value() << " "
        << kappa_.dimensions() << nl
        << "Hardening modulus = " << Hm_.value() << " " << Hm_.dimensions()
        << nl
        << "Initial yield stress = " << Ys0_.value() << " "
        << Ys0_.dimensions() << nl
        << "Linear pressure wave speed = " << pWaveSpeed_.value() << " "
        << pWaveSpeed_.dimensions() << nl
        << "Linear shear wave speed = " << sWaveSpeed_.value() << " "
        << sWaveSpeed_.dimensions() << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

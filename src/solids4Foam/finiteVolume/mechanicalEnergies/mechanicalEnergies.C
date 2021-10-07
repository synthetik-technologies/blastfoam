/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "mechanicalEnergies.H"
#include "fvc.H"

#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mechanicalEnergies, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mechanicalEnergies::mechanicalEnergies
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    externalWork_(0.0),
    externalWorkOldTime_(0.0),
    internalEnergy_(0.0),
    internalEnergyOldTime_(0.0),
    kineticEnergy_(0.0),
    // laplacianSmoothingEnergy_(0.0),
    // laplacianSmoothingEnergyOldTime_(0.0),
    linearBulkViscosityEnergy_(0.0),
    linearBulkViscosityEnergyOldTime_(0.0),
    linearBulkViscosityCoeff_
    (
        dict.lookupOrDefault<scalar>
        (
            "linearBulkViscosityCoeff", 0.06
        )
    ),
    viscousPressurePtr_(),
    epsilonVolPtr_(),
    energiesFilePtr_(),
    curTimeIndex_(-1)
{
    // TODO: read/write energies to allow restart?
    //wip();

    if (Pstream::master())
    {
        Pout<< "Writing energies.dat" << endl;

        energiesFilePtr_.set(new OFstream("energies.dat"));

        energiesFilePtr_()
            << "Time "
            << "External "
            << "Internal "
            << "Kinetic "
            //<< "Smoothing "
            << "Viscosity"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const surfaceScalarField& mechanicalEnergies::viscousPressure
(
    const volScalarField& rho,
    const surfaceScalarField& waveSpeed,
    const volTensorField& gradD
)
{
    if (viscousPressurePtr_.empty())
    {
        viscousPressurePtr_.set
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "viscousPressure",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimPressure, 0.0)
            )
        );
    }

    viscousPressurePtr_() = linearBulkViscosityCoeff_*fvc::interpolate
    (
        rho*fvc::ddt(epsilonVol(gradD))
    )*waveSpeed/mesh_.deltaCoeffs();

    return viscousPressurePtr_();
}


const volScalarField& mechanicalEnergies::epsilonVol
(
    const volTensorField& gradD
)
{
    if (epsilonVolPtr_.empty())
    {
        epsilonVolPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "epsilonVol",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }

    epsilonVolPtr_() = tr(gradD)/3.0;

    return epsilonVolPtr_();
}


void mechanicalEnergies::checkEnergies
(
    const volScalarField& rho,
    const volVectorField& U,
    const volVectorField& D,
    const volVectorField& DD,
    const volSymmTensorField& sigma,
    const volTensorField& gradD,
    const volTensorField& gradDD,
    const surfaceScalarField& waveSpeed,
    const dimensionedVector& g,
    const scalar, // laplacianSmoothCoeff,
    const surfaceScalarField& impKf
)
{
    if (curTimeIndex_ != mesh_.time().timeIndex())
    {
        curTimeIndex_ = mesh_.time().timeIndex();

        // Update old time values
        externalWorkOldTime_ = externalWork_;
        internalEnergyOldTime_ = internalEnergy_;
        //laplacianSmoothingEnergyOldTime_ = laplacianSmoothingEnergy_;
        linearBulkViscosityEnergyOldTime_ = linearBulkViscosityEnergy_;
    }

    // Calculate kinetic energy
    kineticEnergy_ = gSum(0.5*rho.internalField()*mesh_.V()*(U & U));

    // Integrate internal energy using the trapezoidal rule
    internalEnergy_ =
        internalEnergyOldTime_
      + gSum
        (
            DimensionedField<scalar, volMesh>
            (
                mesh_.V()*0.5
               *(
                   sigma.internalField() + sigma.oldTime().internalField()
                ) && symm(gradDD.internalField())
            )
        );

    // Integrate external work energy using the trapezoidal rule
    externalWork_ = externalWorkOldTime_;
    forAll(mesh_.boundary(), patchI)
    {
        if (!mesh_.boundary()[patchI].coupled())
        {
            externalWork_ +=
                gSum
                (
                    (
                        0.5*mesh_.Sf().boundaryField()[patchI]
                      & (
                          sigma.boundaryField()[patchI]
                        + sigma.oldTime().boundaryField()[patchI]
                        )
                    )
                  & DD.boundaryField()[patchI]
                );
        }
    }

    // Include gravity energy
    externalWork_ +=
        gSum
        (
            DimensionedField<scalar, volMesh>
            (
                mesh_.V()*rho.internalField()*g.value() & DD.internalField()
            )
        );

    // Integrate linear bulk viscosity energy using the trapezoidal rule
    if (viscousPressurePtr_.valid())
    {
        linearBulkViscosityEnergy_ =
            linearBulkViscosityEnergyOldTime_
          + gSum
            (
                DimensionedField<scalar, volMesh>
                (
                    fvc::reconstruct
                    (
                        0.5
                       *(
                           viscousPressurePtr_()
                         + viscousPressurePtr_().oldTime()
                        )*mesh_.Sf()
                    )().internalField() && gradDD.internalField()*mesh_.V()
                )
            );
    }

    // Integrate energy dissipated due to Laplacian (Lax-Friedrichs) smoothing
    // term
    //const dimensionedScalar& deltaT = mesh_.time().deltaT();
    //const dimensionedScalar& deltaT0 = mesh_.time().deltaT0();
    // laplacianSmoothingEnergy_ =
    //     laplacianSmoothingEnergyOldTime_
    //   + gSum
    //     (
    //         fvc::reconstruct
    //         (
    //             laplacianSmoothCoeff*0.5*(deltaT + deltaT0)*impKf
    //            *(
    //                fvc::snGrad(U) + fvc::snGrad(U.oldTime())
    //             )*mesh_.magSf()
    //         )().internalField() && gradDD.internalField()*mesh_.V()
    //     );

    // Check the energy imbalance
    // Ideally this should stay less than 1% of the max energy component

    const scalar energyImbalance =
        externalWork_
      - internalEnergy_
      - kineticEnergy_
        //- laplacianSmoothingEnergy_
      - linearBulkViscosityEnergy_;

    const scalar energyImbalancePercent =
        100.0*mag(energyImbalance)/max
        (
            SMALL, max(externalWork_, max(internalEnergy_, kineticEnergy_))
        );

    Info<< "External work = " << externalWork_ << " J" << nl
        << "Internal energy = " << internalEnergy_ << " J" << nl
        << "Kinetic energy = " << kineticEnergy_ << " J" << nl
        //<< "laplacian smoothing energy = "
        //<< laplacianSmoothingEnergy_ << " J"
        << nl
        << "Bulk viscosity energy = " << linearBulkViscosityEnergy_ << " J"
        << nl
        << "Energy imbalance (% of max) = " << energyImbalancePercent << " %"
        << endl;

    if (energyImbalancePercent > 10.0 && debug)
    {
       WarningIn(type() + "::checkEnergies()")
           << "The energy imbalance is greater than 10%" << endl;
       // FatalErrorIn(type() + "::checkEnergies()")
       //     << "The energy imbalance is greater than 10%"
       //     << abort(FatalError);
    }

    // Write energies to file
    if (Pstream::master())
    {
        energiesFilePtr_()
            << mesh_.time().value() << " "
            << externalWork_ << " "
            << internalEnergy_ << " "
            << kineticEnergy_ << " "
            //<< laplacianSmoothingEnergy_ << " "
            << linearBulkViscosityEnergy_
            << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

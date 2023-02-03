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
#include "meshSizeObject.H"

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
    kineticEnergyOldTime_(0.0),
    smoothingEnergy_(0.0),
    smoothingEnergyOldTime_(0.0),
    bulkViscosityEnergy_(0.0),
    bulkViscosityEnergyOldTime_(0.0),
    linearBulkViscosityCoeff_
    (
        "linearBulkViscosityCoeff",
        dimless,
        dict.lookupOrDefault<scalar>
        (
            "linearBulkViscosityCoeff", 0.06
        )
    ),
    quadraticBulkViscosityCoeff_
    (
        "quadraticBulkViscosityCoeff",
        dimless,
        dict.lookupOrDefault<scalar>
        (
            "quadraticBulkViscosityCoeff", 1.2
        )
    ),
    viscousPressurePtr_(),
    energiesFilePtr_(),
    curTimeIndex_(-1)
{
    // TODO: read/write energies to allow restart?
    //wip();

    bool writeEnergies = dict.lookupOrDefault("writeEnergies", false);
    if (Pstream::master() && writeEnergies)
    {
        Pout<< "Writing energies.dat" << endl;

        energiesFilePtr_.set(new OFstream("energies.dat"));

        energiesFilePtr_()
            << "Time "
            << "External "
            << "Internal "
            << "Kinetic "
            << "Smoothing "
            << "Viscosity"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const volScalarField& mechanicalEnergies::viscousPressure
(
    const volScalarField& rho,
    const surfaceScalarField& waveSpeed,
    const volTensorField& gradD
) const
{
    if (viscousPressurePtr_.empty())
    {
        viscousPressurePtr_.set
        (
            new volScalarField
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

    const volScalarField& L(meshSizeObject::New(mesh_).dx(mesh_));
    volScalarField epsilonDot(tr(fvc::ddt(gradD))/3.0);

    viscousPressurePtr_() =
        rho*linearBulkViscosityCoeff_*epsilonDot*fvc::average(waveSpeed)*L;

    epsilonDot.min(0);
    viscousPressurePtr_() +=
        rho*sqr(quadraticBulkViscosityCoeff_*epsilonDot*L);

    return viscousPressurePtr_();
}


const surfaceScalarField& mechanicalEnergies::viscousPressuref
(
    const volScalarField& rho,
    const surfaceScalarField& waveSpeed,
    const volTensorField& gradD
) const
{
    if (viscousPressurefPtr_.empty())
    {
        viscousPressurefPtr_.set
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

    const surfaceScalarField L(1.0/mesh_.deltaCoeffs());
    surfaceScalarField rhof(fvc::interpolate(rho));
    surfaceScalarField epsilonDotf(fvc::interpolate(tr(fvc::ddt(gradD))/3.0));

    viscousPressurefPtr_() =
        rhof*linearBulkViscosityCoeff_*epsilonDotf*waveSpeed*L;

    epsilonDotf.min(0);
    viscousPressurefPtr_() +=
        rhof*sqr(quadraticBulkViscosityCoeff_*epsilonDotf*L);

    return viscousPressurefPtr_();
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
    const momentumStabilisation& stabilisation,
    const dimensionedVector& g
)
{
    Info<<"Energies:" << endl << incrIndent;

    // Store old values
    if (curTimeIndex_ != mesh_.time().timeIndex())
    {
        curTimeIndex_ = mesh_.time().timeIndex();

        // Update old time value
        externalWorkOldTime_ = externalWork_;
        internalEnergyOldTime_ = internalEnergy_;
        kineticEnergyOldTime_ = kineticEnergy_;
        smoothingEnergyOldTime_ = smoothingEnergy_;
        bulkViscosityEnergyOldTime_ = bulkViscosityEnergy_;
    }


    // Write time to output
    if (energiesFilePtr_.valid())
    {
        energiesFilePtr_() << mesh_.time().value();
    }
    scalar energyImbalance = 0;


    // Integrate external work energy using the trapezoidal rule
    {
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
        energyImbalance += externalWork_;

        Info<< indent << "External work = " << externalWork_ << " J" << nl;
        if (energiesFilePtr_.valid())
        {
            energiesFilePtr_()<< token::SPACE << externalWork_;
        }
    }


    // Calculate kinetic energy
    {
        kineticEnergy_ = gSum(0.5*rho.internalField()*mesh_.V()*(U & U));
        energyImbalance -= kineticEnergy_;

        Info<< indent << "Kinetic energy = " << kineticEnergy_ << " J" << nl;
        if (energiesFilePtr_.valid())
        {
            energiesFilePtr_()<< " " << kineticEnergy_;
        }
    }

    // Integrate internal energy using the trapezoidal rule
    {
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
        energyImbalance -= internalEnergy_;

        Info<< indent << "Internal energy = " << internalEnergy_ << " J" << nl;
        if (energiesFilePtr_.valid())
        {
            energiesFilePtr_()<< token::SPACE << internalEnergy_;
        }
    }

    // Integrate linear bulk viscosity energy using the trapezoidal rule
    if (viscousPressurefPtr_.valid() || viscousPressurePtr_.valid())
    {
        if (viscousPressurePtr_.valid())
        {
            bulkViscosityEnergy_ =
                bulkViscosityEnergyOldTime_
              + gSum
                (
                    DimensionedField<scalar, volMesh>
                    (
                        // fvc::grad
                        (
                            0.5
                           *(
                                viscousPressurePtr_()
                              + viscousPressurePtr_().oldTime()
                            )*tensor::I
                        )().internalField() && (gradDD.internalField()*mesh_.V())
                    )
                );
        }
        else if (viscousPressurefPtr_.valid())
        {
            bulkViscosityEnergy_ =
                bulkViscosityEnergyOldTime_
              + gSum
                (
                    DimensionedField<scalar, volMesh>
                    (
                        fvc::reconstruct
                        (
                            0.5
                           *(
                                viscousPressurefPtr_()
                              + viscousPressurefPtr_().oldTime()
                            )*mesh_.Sf()
                        )().internalField() && (gradDD.internalField()*mesh_.V())
                    )
                );
        }
        energyImbalance -= bulkViscosityEnergy_;

        Info<< indent << "Bulk viscosity energy = "
            << bulkViscosityEnergy_ << " J" << nl;
        if (energiesFilePtr_.valid())
        {
            energiesFilePtr_()<< token::SPACE << bulkViscosityEnergy_;
        }
    }

    // Integrate energy dissipated due to Laplacian (Lax-Friedrichs) smoothing
    // term
    if (stabilisation.inUse())
    {
        const dimensionedScalar& deltaT01 =
            0.5*(mesh_.time().deltaT() + mesh_.time().deltaT0());
        smoothingEnergy_ =
            smoothingEnergyOldTime_
          + stabilisation.energy
            (
                U,
                (mesh_.lookupObject<surfaceScalarField>("impKf")*deltaT01)(),
                gradDD
            );
        energyImbalance -= smoothingEnergy_;

        Info<< indent << "Smoothing energy = "
            << smoothingEnergy_ << " J" << nl;

        if (energiesFilePtr_.valid())
        {
            energiesFilePtr_()<< token::SPACE << smoothingEnergy_;
        }
    }

    // Check the energy imbalance
    // Ideally this should stay less than 1% of the max energy component

    const scalar energyImbalancePercent =
        100.0*mag(energyImbalance)/max
        (
            SMALL, max(externalWork_, max(internalEnergy_, kineticEnergy_))
        );

    Info<< indent << "Energy imbalance (% of max) = "
        << energyImbalancePercent << " %" << endl;

    if (energyImbalancePercent > 10.0 && debug)
    {
       WarningInFunction
           << "The energy imbalance is greater than 10%" << endl;
    }

    // Write energies to file
    if (energiesFilePtr_.valid())
    {
        energiesFilePtr_()<< endl;
    }

    Info<< endl << decrIndent;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

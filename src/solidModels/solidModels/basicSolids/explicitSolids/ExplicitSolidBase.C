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

#include "ExplicitSolidBase.H"
#include "solidTractionFvPatchVectorField.H"
#include "wedgeFvPatch.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * * //

template<class IncrementalSolid>
void Foam::solidModels::ExplicitSolidBase<IncrementalSolid>::relax()
{
    // Re-read the dictionary
    relaxation_.read
    (
        this->solidModelDict().optionalSubDict("relaxation")
    );

    //- Relax, if wanted
    relaxation_.relax(this->U(), this->DD(), a_, this->rho());
}


template<class IncrementalSolid>
void Foam::solidModels::ExplicitSolidBase<IncrementalSolid>::correctUBCs
(
    volVectorField& U
)
{
    U.correctBoundaryConditions();
    volVectorField::Boundary& bU = U.boundaryFieldRef();
    const volVectorField::Boundary& bDD = this->DD().boundaryField();
    forAll(bU, patchi)
    {
        if (isA<solidTractionFvPatchVectorField>(bDD[patchi]))
        {
            const fvPatch& patch = this->mesh().boundary()[patchi];

            fvPatchVectorField& pU = bU[patchi];
            const solidTractionFvPatchVectorField& pDD =
                dynamicCast<const solidTractionFvPatchVectorField>(bDD[patchi]);
            vectorField n(this->nf(patch));
            tensorField nn(n*n);

            tensorField St
            (
                nn/this->wavespeed_.boundaryField()[patchi]
            + (I - nn)/this->sWavespeed_.boundaryField()[patchi]
            );

            pU =
                pU.internalField()
              + (
                    St
                  & (
                        pDD.traction() - pDD.pressure()*n
                      - this->Pn(patch)
                    )
                )/this->rho().boundaryField()[patchi];
        }
        else if (bDD[patchi].fixesValue())
        {
            bU[patchi] = bDD[patchi]/this->mesh().time().deltaTValue();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IncrementalSolid>
Foam::solidModels::ExplicitSolidBase<IncrementalSolid>::ExplicitSolidBase
(
    const word& type,
    dynamicFvMesh& mesh,
    const nonLinearGeometry::nonLinearType nonLinear,
    const bool isSolid
)
:
    IncrementalSolid(type, mesh, nonLinear, isSolid),
    wavespeed_
    (
        IOobject
        (
            "wavespeed",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimVelocity, Zero)
    ),
    sWavespeed_
    (
        IOobject
        (
            "sWavespeed",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimVelocity, Zero)
    ),
    energies_(mesh, this->solidModelDict()),
    a_
    (
        IOobject
        (
            "a",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity/dimTime, Zero),
        extrapolatedCalculatedFvPatchScalarField::typeName
    ),
    relaxation_(this->solidModelDict().optionalSubDict("relaxation"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IncrementalSolid>
void Foam::solidModels::ExplicitSolidBase<IncrementalSolid>::solveMomentum()
{
    Info<< "Solving the momentum equation" << endl;

    this->enforceLinear() = false;

    do
    {
        // Central difference scheme
        const dimensionedScalar& deltaT = this->time().deltaT();
        const dimensionedScalar deltaT01(0.5*(deltaT + this->time().deltaT0()));

        // Compute the velocity
        // Note: this is the velocity at the middle of the time-step
        this->U() = this->U().oldTime() + deltaT01*a_.oldTime();

        // Compute change in displacement
        this->DD().ref() = deltaT*this->U()();

        // Enforce any cell displacements
        if (this->setCellDisps().cellIDs().size())
        {
            vectorField& DDI = this->DD();
            vectorField& UI = this->U();
            const vectorField& DOld = this->D().oldTime();

            const labelList& cells = this->setCellDisps().cellIDs();
            const vectorField& cellDs = this->setCellDisps().cellDisps();

            forAll(cells, i)
            {
                const label celli = cells[i];
                DDI[celli] = cellDs[i] - DOld[celli];
                UI[celli] = DDI[celli]/deltaT.value();
            }
        }
        this->DD().correctBoundaryConditions();
        this->U().boundaryFieldRef() =
            this->DD().boundaryField()/deltaT.value();
        this->correctUBCs(this->U());

        relax();

        // Update displacement
        this->D() == this->D().oldTime() + this->U()*deltaT;

        // Update the stress field based on the latest D field
        this->update();

        // Compute acceleration
        // Note the inclusion of a linear bulk viscosity pressure term to
        // dissipate high frequency energies, and a Rhie-Chow term to
        // avoid checker-boarding
        tmp<volVectorField> stab
        (
            this->stabilisation().stabilisation
            (
                this->U(),
                fvc::grad(this->U())(),
               (deltaT01*this->impKf_)()
            )
        );

        surfaceVectorField tractionSf(this->tractionSf());
        // surfaceVectorField::Boundary& btractionSf = tractionSf.boundaryFieldRef();
        // const volVectorField::Boundary& bD = this->solutionD().boundaryField();
        // forAll(btractionSf, patchi)
        // {
        //     if (isA<solidTractionFvPatchVectorField>(bD[patchi]))
        //     {
        //         const fvPatch& patch = this->mesh().boundary()[patchi];
        //         const solidTractionFvPatchVectorField& pD =
        //             dynamicCast<const solidTractionFvPatchVectorField>(bD[patchi]);
        //
        //         // Set boundary traction
        //         btractionSf[patchi] =
        //             (
        //                 pD.traction() - this->nf(patch)*pD.pressure()
        //               // + (this->nf(patch) & this->sigma(patch))
        //             )*patch.magSf();
        //     }
        // }
        a_ =
            (
                fvc::div(tractionSf)
                // this->divStress()
              + fvc::div
                (
                    this->mesh().Sf()*energies_.viscousPressure
                    (
                        this->rho(),
                        wavespeed_,
                        this->gradD()
                    )
                )

                // This corresponds to Lax–Friedrichs smoothing
              + stab()
            )/this->rho()
          + this->g();
        a_.correctBoundaryConditions();

        // Check energies
        energies_.checkEnergies
        (
            this->rho(),
            this->U(),
            this->D(),
            this->DD(),
            this->sigma(),
            this->gradD(),
            this->gradDD(),
            stab(),
            this->g()
        );

    } while (this->mesh().update());

    // Mesh update loop
    // do
    // {
    //     // Central difference scheme
    //     const dimensionedScalar& deltaT = this->time().deltaT();
    //     const dimensionedScalar deltaT01(0.5*(deltaT + this->time().deltaT0()));
    //
    //     // Compute the velocity
    //     // Note: this is the velocity at the middle of the time-step
    //     this->U() = this->U().oldTime() + deltaT01*a_.oldTime();
    //
    //     // Compute change in displacement
    //     this->DD() = deltaT*this->U();
    //
    //     // Enforce any cell displacements
    //     if (this->setCellDisps().cellIDs().size())
    //     {
    //         vectorField& DDI = this->DD();
    //         vectorField& UI = this->U();
    //         const vectorField& DOld = this->D().oldTime();
    //
    //         const labelList& cells = this->setCellDisps().cellIDs();
    //         const vectorField& cellDs = this->setCellDisps().cellDisps();
    //
    //         forAll(cells, i)
    //         {
    //             const label celli = cells[i];
    //             DDI[celli] = cellDs[i] - DOld[celli];
    //             UI[celli] = DDI[celli]/deltaT.value();
    //         }
    //     }
    //     this->DD().correctBoundaryConditions();
    //     this->U().boundaryFieldRef() =
    //         this->DD().boundaryField()/this->mesh().time().deltaTValue();
    //     // this->correctUBCs(this->U());
    //
    //     relax();
    //
    //     // // Compute change in displacement
    //     // this->DD() == deltaT*this->U();
    //
    //     // Update displacement
    //     this->D() == this->D().oldTime() + this->DD();//*deltaT;
    //
    //     // Update the stress field based on the latest D field
    //     this->update();
    //
    //     surfaceVectorField tractionSf(this->tractionSf());
    //     surfaceVectorField::Boundary& btractionSf = tractionSf.boundaryFieldRef();
    //     const volVectorField::Boundary& bD = this->solutionD().boundaryField();
    //     forAll(btractionSf, patchi)
    //     {
    //         if (isA<solidTractionFvPatchVectorField>(bD[patchi]))
    //         {
    //             const fvPatch& patch = this->mesh().boundary()[patchi];
    //             const solidTractionFvPatchVectorField& pD =
    //                 dynamicCast<const solidTractionFvPatchVectorField>(bD[patchi]);
    //
    //             // Set boundary traction
    //             btractionSf[patchi] =
    //                 (
    //                     pD.traction() - this->nf(patch)*pD.pressure()
    //                   + (this->nf(patch) & this->sigma(patch))
    //                 )*patch.magSf();
    //         }
    //     }
    //
    //     // Compute acceleration
    //     // Note the inclusion of a linear bulk viscosity pressure term to
    //     // dissipate high frequency energies, and a Rhie-Chow term to
    //     // avoid checker-boarding
    //     tmp<volVectorField> stab
    //     (
    //         this->stabilisation().stabilisation
    //         (
    //             this->U(),
    //             fvc::grad(this->U())(),
    //            (deltaT01*this->impKf_)()
    //         )
    //     );
    //
    //     a_ =
    //         (
    //             fvc::div(tractionSf)
    //           + fvc::div
    //             (
    //                 this->mesh().Sf()*energies_.viscousPressure
    //                 (
    //                     this->rho(),
    //                     wavespeed_,
    //                     this->gradD()
    //                 )
    //             )
    //
    //             // This corresponds to Lax–Friedrichs smoothing
    //           + stab()
    //         )/this->rho()
    //       + this->g();
    //     a_.correctBoundaryConditions();
    //
    //     // Check energies
    //     energies_.checkEnergies
    //     (
    //         this->rho(),
    //         this->U(),
    //         this->D(),
    //         this->DD(),
    //         this->sigma(),
    //         this->gradD(),
    //         this->gradDD(),
    //         stab(),
    //         this->g()
    //     );
    //
    // } while (this->mesh().update());
}


template<class IncrementalSolid>
Foam::scalar
Foam::solidModels::ExplicitSolidBase<IncrementalSolid>::CoNum() const
{
    // waveSpeed = cellWidth/deltaT
    // So, deltaT = cellWidth/waveVelocity == (1.0/deltaCoeff)/waveSpeed
    // In the current discretisation, information can move two cells per
    // time-step. This means that we use 1/(2*d) == 0.5*deltaCoeff when
    // calculating the required stable time-step
    // i.e.e deltaT = (1.0/(0.5*deltaCoeff)/waveSpeed
    // For safety, we should use a time-step smaller than this e.g. Abaqus uses
    // 1/sqrt(2)*stableTimeStep: we will default to this value

    const fvMesh& mesh = this->mesh();
    surfaceScalarField amaxSf(this->wavespeed()*mesh.magSf());

    // Remove wave speed from wedge boundaries
    forAll(amaxSf.boundaryField(), patchi)
    {
        if (isA<wedgeFvPatch>(mesh.boundary()[patchi]))
        {
            amaxSf.boundaryFieldRef() = Zero;
        }
    }
    scalarField sumAmaxSf
    (
        fvc::surfaceSum(amaxSf)().primitiveField()
    );
    return
        0.5*gMax(sumAmaxSf/mesh.V().field())*mesh.time().deltaTValue();
}


template<class IncrementalSolid>
Foam::scalar
Foam::solidModels::ExplicitSolidBase<IncrementalSolid>::maxCoNum() const
{
    return
        this->mesh().time().controlDict().lookupOrDefault
        (
            "maxCo",
            scalar(0.7071)
        );
}

// ************************************************************************* //

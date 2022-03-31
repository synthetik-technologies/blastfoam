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

#include "blastFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "dynamicBlastFvMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "wedgeFvPatch.H"
#include "fluidSolidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(blastFluid, 0);
addToRunTimeSelectionTable(physicsModel, blastFluid, fluid);
addToRunTimeSelectionTable(fluidModel, blastFluid, dictionary);


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

blastFluid::blastFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region, true)
#ifdef OPENFOAMFOUNDATION
    , integrator_(timeIntegrator::New(mesh())),
    fluid_(compressibleSystem::New(mesh())),
    p_(fluid_->p()),
    T_(fluid_->T()),
    curTimeIndex_(-1)
{
    integrator_->addSystem(fluid_());

    fluid_->update();
    integrator_->clear();
}
#else
{}
#endif

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> blastFluid::patchViscousForce
(
    const label patchID
) const
{
#ifdef OPENFOAMFOUNDATION
    if (!fluid_->inviscid())
    {
        return
            mesh().boundary()[patchID].Sf()
          & fluid_->turbulence().devTau()().boundaryField()[patchID];
    }
#endif
    return tmp<vectorField>
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

}


tmp<scalarField> blastFluid::patchPressureForce
(
    const label patchID
) const
{
    // Pressure here is already in Pa
#ifdef OPENFOAMFOUNDATION
    return tmp<scalarField>
    (
        new scalarField(p_.boundaryField()[patchID])
    );
#else
    return tmp<scalarField>
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );
#endif

}


bool blastFluid::evolve()
{
#ifdef OPENFOAMFOUNDATION
    Info<< "Evolving fluid model: " << this->type() <<nl<< endl;

    dynamicFvMesh& mesh = this->mesh();

    if (fluidModel::fsiMeshUpdate())
    {
        // The FSI interface is in charge of calling mesh.update()
        fluidModel::fsiMeshUpdateChanged();
    }
    else
    {
        mesh.update();
    }

    //- ODE integration
    integrator_->integrate();

    // Info variable values
    Info<< "max(p): " << max(p_).value()
        << ", min(p): " << min(p_).value() << endl;
    Info<< "max(T): " << max(T_).value()
        << ", min(T): " << min(T_).value() << endl;

    //- Clear tmp fields
    integrator_->clear();
#else
    NotImplemented;
#endif
    return false;
}


void blastFluid::setDeltaT(Time& runTime)
{
    // Refine mesh before setting the time step
    // also removes problems with incorrectly sized global patches
    refineMesh(mesh());

    if (!runTime.controlDict().lookupOrDefault("adjustTimeStep", true))
    {
        return;
    }

    // Calculate the maximum Courant number
    surfaceScalarField amaxSf
    (
        fvc::interpolate(fluid_->speedOfSound())*mesh().magSf()
    );
    // Remove wave speed from wedge boundaries
    forAll(amaxSf.boundaryField(), patchi)
    {
        if (isA<wedgeFvPatch>(mesh().boundary()[patchi]))
        {
            amaxSf.boundaryFieldRef()[patchi] = Zero;
        }
    }
    amaxSf += mag(fvc::flux(fluid_->U()));

    scalarField sumAmaxSf
    (
        fvc::surfaceSum(amaxSf)().primitiveField()
    );

    scalar CoNum =
        0.5
       *gMax(sumAmaxSf/mesh().V().field())
       *mesh().time().deltaTValue();

    scalar maxCo
    (
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.5)
    );
    scalar maxDeltaT
    (
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT)
    );

    scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
    scalar deltaTFact =
        min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaT().value(),
            maxDeltaT
        )
    );

    Info<< "Max Courant Number = " << CoNum << endl;
    Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "phaseCompressibleSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseCompressibleSystem, 0);
    defineRunTimeSelectionTable(phaseCompressibleSystem, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseCompressibleSystem::phaseCompressibleSystem
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 0.0)
    ),
    U_
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    rhoU_
    (
        IOobject
        (
            "rhoU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho_*U_
    ),
    rhoE_
    (
        IOobject
        (
            "rhoE",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho_*(e_ + 0.5*magSqr(U_))
    ),
    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh
        ),
        fvc::flux(U_)
    ),
    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity*dimVelocity*dimArea, 0.0)
    ),
    rhoUPhi_
    (
        IOobject
        (
            "rhoUPhi",
            mesh.time().timeName(),
            mesh
        ),
        fvc::interpolate(rhoU_)*phi_
    ),
    rhoEPhi_
    (
        IOobject
        (
            "rhoEPhi",
            mesh.time().timeName(),
            mesh
        ),
        fvc::interpolate(rhoE_)*phi_
    ),
    fluxScheme_(fluxScheme::New(mesh))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseCompressibleSystem::~phaseCompressibleSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseCompressibleSystem::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    if (oldIs_[stepi - 1] != -1)
    {
        rhoUOld_[oldIs_[stepi - 1]] = rhoU_;
        rhoEOld_[oldIs_[stepi - 1]] = rhoE_;
    }
    volVectorField rhoUOld(ai[stepi - 1]*rhoU_);
    volScalarField rhoEOld(ai[stepi - 1]*rhoE_);

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = oldIs_[i];
        if (fi != -1 && ai[fi] != 0)
        {
            rhoUOld += ai[fi]*rhoUOld_[fi];
            rhoEOld += ai[fi]*rhoEOld_[fi];
        }
    }

    volVectorField deltaRhoU(fvc::div(rhoUPhi_));
    volScalarField deltaRhoE(fvc::div(rhoEPhi_));
    if (deltaIs_[stepi - 1] != -1)
    {
        deltaRhoU_[deltaIs_[stepi - 1]] = deltaRhoU;
        deltaRhoE_[deltaIs_[stepi - 1]] = deltaRhoE;
    }
    deltaRhoU *= bi[stepi - 1];
    deltaRhoE *= bi[stepi - 1];

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            deltaRhoU += bi[fi]*deltaRhoU_[fi];
            deltaRhoE += bi[fi]*deltaRhoE_[fi];
        }
    }


    dimensionedScalar dT = rho_.time().deltaT();
    vector solutionDs((vector(rho_.mesh().solutionD()) + vector::one)/2.0);

    rhoU_ = cmptMultiply(rhoUOld - dT*deltaRhoU, solutionDs);
    rhoU_.correctBoundaryConditions();

    rhoE_ = rhoEOld - dT*deltaRhoE;
    rhoE_.correctBoundaryConditions();
}


void Foam::phaseCompressibleSystem::setODEFields
(
    const label nSteps,
    const boolList& storeFields,
    const boolList& storeDeltas
)
{
    oldIs_.resize(nSteps);
    deltaIs_.resize(nSteps);
    label fi = 0;
    for (label i = 0; i < nSteps; i++)
    {
        if (storeFields[i])
        {
            oldIs_[i] = fi;
            rhoUOld_.append
            (
                new volVectorField
                (
                    IOobject
                    (
                        IOobject::groupName(rhoU_.name(), Foam::name(i)),
                        rho_.time().timeName(),
                        rho_.mesh()
                    ),
                    rhoU_
                )
            );
            rhoEOld_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(rhoE_.name(), Foam::name(i)),
                        rho_.time().timeName(),
                        rho_.mesh()
                    ),
                    rhoE_
                )
            );

            fi++;
        }
        else
        {
            oldIs_[i] = -1;
        }
    }
    fi = 0;
    for (label i = 0; i < nSteps; i++)
    {
        if (storeDeltas[i])
        {
            deltaIs_[i] = fi;
            deltaRhoU_.append
            (
                new volVectorField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            rhoU_.name() + "Delta", Foam::name(i)
                        ),
                        rho_.time().timeName(),
                        rho_.mesh()
                    ),
                    rho_.mesh(),
                    dimensionedVector("0", rhoU_.dimensions()/dimTime, Zero)
                )
            );
            deltaRhoE_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            rhoE_.name() + "Delta", Foam::name(i)
                        ),
                        rho_.time().timeName(),
                        rho_.mesh()
                    ),
                    rho_.mesh(),
                    dimensionedScalar("0", rhoE_.dimensions()/dimTime, 0.0)
                )
            );
            fi++;
        }
        else
        {
            deltaIs_[i] = -1;
        }
    }
}

// ************************************************************************* //

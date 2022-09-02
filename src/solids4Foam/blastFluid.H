/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
31-03-2022 Synthetik Appplied Technologies: Added blastFoam based fluidModel
-------------------------------------------------------------------------------

License
    This file is a derived work of OpenFOAM.

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

Class
    blastFluid

Description
    Transient solver for trans-sonic/supersonic, laminar flow of a
    compressible liquid. Assumes a linear equation of state between
	density and pressure related by the bulk modulus.

Author
    Iago Lessa de Oliveira, FEIS/UNESP/UCD.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    blastFluid.C

\*---------------------------------------------------------------------------*/

#ifndef blastFluid_H
#define blastFluid_H

#include "fluidModel.H"
#include "compressibleSystem.H"
#include "timeIntegrator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class blastFluid Declaration
\*---------------------------------------------------------------------------*/

class blastFluid
:
    public fluidModel
{
#ifdef OPENFOAMFOUNDATION
    // Private data

        //- Time integrator
        autoPtr<timeIntegrator> integrator_;

        // Compressible system
        autoPtr<compressibleSystem> fluid_;

        //- Pressure
        const volScalarField& p_;

        //- Temperature
        const volScalarField& T_;

        //- Current time index
        label curTimeIndex_;

#endif
        //- Disallow default bitwise copy construct
        blastFluid(const blastFluid&);

        //- Disallow default bitwise assignment
        void operator=(const blastFluid&);


public:

    //- Runtime type information
    TypeName("blast");

    // Constructors

        //- Construct from components
        blastFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~blastFluid()
        {}


    // Member Functions

        // Access

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> patchViscousForce
            (
                const label patchID
            ) const;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> patchPressureForce
            (
                const label patchID
            ) const;

            //- Return const velocity field
            virtual const volVectorField& U() const
            {
                return fluid_->U();
            }

            //- Return velocity field
            virtual volVectorField& U()
            {
                NotImplemented;
                return fluidModel::U();
            }

            //- Return const pressure field
            virtual const volScalarField& p() const
            {
                return p_;
            }

            //- Return pressure field
            virtual volScalarField& p()
            {
                return fluid_->thermo().p();
            }

            //- Return const velocity gradient field
            virtual const volTensorField& gradU() const
            {
                NotImplemented;
                return fluidModel::gradU();
            }

            //- Return velocity gradient field
            virtual volTensorField& gradU()
            {
                NotImplemented;
                return fluidModel::gradU();
            }

            //- Return const pressure gradient field
            virtual const volVectorField& gradp() const
            {
                NotImplemented;
                return fluidModel::gradp();
            }

            //- Return pressure gradient field
            virtual volVectorField& gradp()
            {
                NotImplemented;
                return fluidModel::gradp();
            }

            //- Return const reference to flux field
            virtual const surfaceScalarField& phi() const
            {
                return fluid_->phi();
            }

            //- Return reference to flux field
            virtual surfaceScalarField& phi()
            {
                NotImplemented;
                return fluidModel::phi();
            }

        // Edit

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime);

            //- Evolve the fluid model
            virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

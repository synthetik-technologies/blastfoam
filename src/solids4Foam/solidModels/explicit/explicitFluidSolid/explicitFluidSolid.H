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

Class
    explicitFluidSolid

Description
    Version of the linGeomTotalDispSolid solid model, where a fully explicit
    solution algorithm is used.

    The explicit implementation uses central differencing, where a linear bulk
    viscosity pressure term dissipates high frequencies.

    Details of the central differencing solution algorithm can be found in, e.g.
    - Belytschko, Liu, Moran, 2000, Box 6.1
    - Ding and Ye, Laser shock peening, Performance and process simulation,
      2006, Section 3.5.2

    The linear bulk viscosity term is described in the Abaqus documentation.

    A Jameson-Schmidt-Turkel (JST) 4th order diffusion term is used for
    stabilisation.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    explicitFluidSolid.C

\*---------------------------------------------------------------------------*/

#ifndef explicitFluidSolid_H
#define explicitFluidSolid_H

#include "solidModel.H"
#include "updatedLagSolid.H"
#include "incrementalSolid.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "mechanicalEnergies.H"

#include "timeIntegrationSystem.H"
#include "fluidBlastThermo.H"
#include "fluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
            Class explicitFluidSolid Declaration
\*---------------------------------------------------------------------------*/

class explicitFluidSolid
:
    public updatedLagSolid<incrementalSolid>,
    public timeIntegrationSystem
{
    // Private data

        //- Reference to the time integrator
        timeIntegrator& integrator_;

        //- ThermodynamicModel
        fluidBlastThermo& thermo_;

        //- Total mass
        volScalarField& rho_;

        //- Pressure
        volScalarField& p_;

        //- Temperature
        volScalarField& T_;

        //- Internal energy
        volScalarField& e_;

        //- Velocity
        volVectorField& U_;

    //- Conserved variables

        //- Momentum
        volVectorField rhoU_;

        //- Energy
        volScalarField rhoE_;


    //- Fluxes

        //- Volumetric flux
        surfaceScalarField phi_;

        //- Total mass flux
        surfaceScalarField rhoPhi_;

        //- Momentum flux
        surfaceVectorField rhoUPhi_;

        //- Energy flux
        surfaceScalarField rhoEPhi_;

        //- Flux scheme
        autoPtr<fluxScheme> fluxScheme_;

        // Scale factor for Lax-Lax–Friedrichs smoothing term
        const scalar LFScaleFactor_;

        // Scale factor for Jameson-Schmidt-Turkel smoothing term
        const scalar JSTScaleFactor_;

        //- Energy calculation
        mechanicalEnergies energies_;


    // Private Member Functions

        //- Decode primative variables
        virtual void decode();

        //- Encode conserved variables
        virtual void encode();

        //- Disallow default bitwise copy construct
        explicitFluidSolid
        (
            const explicitFluidSolid&
        );

        //- Disallow default bitwise assignment
        void operator=(const explicitFluidSolid&);

public:

    //- Runtime type information
    TypeName("explicitFluid");

    // Constructors

        //- Construct from components
        explicitFluidSolid(dynamicFvMesh& mesh);

    // Destructor

        virtual ~explicitFluidSolid()
        {}


    // Member Functions

        // Edit

            //- Update fluxes
            virtual void update();

            //- Solve sub-step stepi
            virtual void solve();

            //- Post update
            virtual void postUpdate();

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime);

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Return the current Courant number
            virtual scalar CoNum() const;

            //- Return the max Courant number
            virtual scalar maxCoNum() const;

            //- End of the program
            virtual void end()
            {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

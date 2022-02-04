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
    explicitNonLinGeomTotalLagSolid

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
    explicitNonLinGeomTotalLagSolid.C

\*---------------------------------------------------------------------------*/

#ifndef explicitNonLinGeomTotalLagSolid_H
#define explicitNonLinGeomTotalLagSolid_H

#include "solidModel.H"
#include "totalLagSolid.H"
#include "totalDispSolid.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "mechanicalEnergies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
            Class explicitNonLinGeomTotalLagSolid Declaration
\*---------------------------------------------------------------------------*/

class explicitNonLinGeomTotalLagSolid
:
    public totalLagSolid<totalDispSolid>
{
    // Private data

        //- Scale factor for Rhie-Chow smoothing term
        //const scalar RhieChowScaleFactor_;

        // Scale factor for Lax-Lax–Friedrichs smoothing term
        const scalar LFScaleFactor_;

        // Scale factor for Jameson-Schmidt-Turkel smoothing term
        const scalar JSTScaleFactor_;

        //- Wave speed field
        surfaceScalarField waveSpeed_;

        //- Energy calculation
        mechanicalEnergies energies_;

        //- Acceleration
        volVectorField a_;


    // Private Member Functions

        //- Update the stress field
        void updateStress();

        //- Smooth the hydrostatic pressure field
        void smoothPressure();

        //- Return smoothed U old time based on Lax-Friedlichs method
        tmp<volVectorField> smoothUOldTime();

        //- Disallow default bitwise copy construct
        explicitNonLinGeomTotalLagSolid
        (
            const explicitNonLinGeomTotalLagSolid&
        );

        //- Disallow default bitwise assignment
        void operator=(const explicitNonLinGeomTotalLagSolid&);

public:

    //- Runtime type information
    TypeName("explicitNonLinearTotalLagrangian");

    // Constructors

        //- Construct from components
        explicitNonLinGeomTotalLagSolid(dynamicFvMesh& mesh);

    // Destructor

        virtual ~explicitNonLinGeomTotalLagSolid()
        {}


    // Member Functions

        // Edit

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

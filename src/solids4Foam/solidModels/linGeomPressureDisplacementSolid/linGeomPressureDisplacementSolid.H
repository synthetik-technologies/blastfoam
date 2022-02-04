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
    linGeomPressureVelocitySolid

Description
    Mathematical model where linear geometry is assumed i.e. small strains and
    small rotations are assumed, and the total displacement and hydrostatic
    pressure are the primary unknowns.

    Rhie-Chow stabilisation is applied directly to the pressure field.

    The stress is calculated by the run-time selectable mechanical law.

Author
    Philip Cardiff, UCD.  All rights reserved.
    Andrew Whelan, UCD.  All rights reserved.

SourceFiles
    linGeomPressureVelocitySolid.C

\*---------------------------------------------------------------------------*/

#ifndef linGeomPressureDisplacementSolid_H
#define linGeomPressureDisplacementSolid_H

#include "linSolid.H"
#include "totalDispSolid.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class linGeomTotalDispSolid Declaration
\*---------------------------------------------------------------------------*/

class linGeomPressureDisplacementSolid
:
    public linSolid<totalDispSolid>
{
    // Private data

        //- Bulk modulus field
        const volScalarField k_;

        //- Reciprocal of bulk modulus field, stored for efficiency
        const volScalarField rK_;

        //- Hydrostatic pressure
        volScalarField p_;

        //- Gradient of hydrostatic pressure
        volVectorField gradp_;

        //- Solve pressure equation implicitly
        const Switch solvePressureEquationImplicitly_;

        //- Scale factor for Rhie-Chow smoothing term in pressure calculation
        const scalar pressureRhieChowScaleFac_;

        //- Number of correctors for inner pressure loop
        const int nInnerCorr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        linGeomPressureDisplacementSolid
        (
            const linGeomPressureDisplacementSolid&
        );

        //- Disallow default bitwise assignment
        void operator=(const linGeomPressureDisplacementSolid&);


public:

    //- Runtime type information
    TypeName("linearPressureDisplacement");

    // Constructors

        //- Construct from components
        linGeomPressureDisplacementSolid(dynamicFvMesh& mesh);

    // Destructor

        virtual ~linGeomPressureDisplacementSolid()
        {}


    // Member Functions

        //- Evolve the solid solver and solve the mathematical model
        virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

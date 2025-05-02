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
    weakThermalLinearSolid

Description
    Mathematical model where linear geometry is assumed i.e. small strains and
    small rotations are assumed, and the heat equation is solved.

    The stress is calculated by the run-time selectable mechanical law.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    weakThermalLinearSolid.C

\*---------------------------------------------------------------------------*/

#ifndef weakThermalLinearSolid_H
#define weakThermalLinearSolid_H

#include "linearTotalDisplacementSolid.H"
#include "thermalModel.H"
#include "SolverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class weakThermalLinearSolid Declaration
\*---------------------------------------------------------------------------*/

class weakThermalLinearSolid
:
    public linearTotalDisplacementSolid
{
    // Private data

        //- The product of density and specific heat
        volScalarField rhoC_;

        //- Thermal conductivity
        volScalarField k_;

        //- Temperature
        volScalarField& T_;

        //- Temperature gradient
        volVectorField gradT_;

        //- Absolute temperature tolerance (in degrees)
        scalar absTTol_;

    // Private Member Functions

        //- Check convergence
        bool converged
        (
            const int iCorr,
            const SolverPerformance<scalar>& solverPerfT,
            const volScalarField& T
        );

    //- Disallow default bitwise copy construct
        weakThermalLinearSolid(const weakThermalLinearSolid&);

        //- Disallow default bitwise assignment
        void operator=(const weakThermalLinearSolid&);

        virtual void readDict();


public:

    //- Runtime type information
    TypeName("weakThermalLinear");

    // Constructors

        //- Construct from components
        weakThermalLinearSolid(fvMesh& mesh);


    // Destructor

        virtual ~weakThermalLinearSolid()
        {}


    // Member Functions

        // Access

            //- Return const access to the cell temperature field
            const volScalarField& T() const
            {
                return T_;
            }

            //- Return non-const access to the cell temperature field
            volScalarField& T()
            {
                return T_;
            }

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Write fields
            virtual bool write(const bool write = true) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

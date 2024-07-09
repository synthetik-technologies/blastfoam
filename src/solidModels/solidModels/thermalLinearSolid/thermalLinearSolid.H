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
    thermalLinearSolid

Description
    Mathematical model where linear geometry is assumed i.e. small strains and
    small rotations are assumed.

    The heat equation and momentum equation are solved in a strongly coupled
    fashion i.e. just like in solidDisplacementFoam.

    The stress is calculated by the run-time selectable mechanical law.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    thermalLinearSolid.C

\*---------------------------------------------------------------------------*/

#ifndef thermalLinearSolid_H
#define thermalLinearSolid_H

#include "LinearGeomSolid.H"
#include "totalDisplacementSolids.H"
#include "thermalModel.H"
#include "SolverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class thermalLinearSolid Declaration
\*---------------------------------------------------------------------------*/

class thermalLinearSolid
:
    public LinearGeomSolid<totalDisplacementSolid>
{
    // Private data

        //- Absolute energy tolerance
        scalar absHETol_;


    // Private Member Functions

        //- Check if the momentum equation is converged
        bool converged
        (
            const int iCorr,
            const SolverPerformance<vector>& solverPerfD,
            const SolverPerformance<scalar>& solverPerfHE,
            const volVectorField& D,
            const volScalarField& he
        );

        //- Disallow default bitwise copy construct
        thermalLinearSolid(const thermalLinearSolid&);

        //- Disallow default bitwise assignment
        void operator=(const thermalLinearSolid&);

        virtual void readDict();


public:

    //- Runtime type information
    TypeName("thermalLinear");

    // Constructors

        //- Construct from components
        thermalLinearSolid(dynamicFvMesh& mesh);


    // Destructor

        virtual ~thermalLinearSolid()
        {}


    // Member Functions

        // Edit

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

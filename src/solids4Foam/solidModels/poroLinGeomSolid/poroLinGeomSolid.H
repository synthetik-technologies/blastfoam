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
    poroLinGeomSolid

Description
    Mathematical model for poro-elasticity, where linear geometry is assumed
    i.e. small strains and small rotations are assumed.

    The poro-pressure equation is solved.

    The stress is calculated by the run-time selectable mechanical law.

    The approach is based on the elastoPlasticBiotFoam in the miniGeotechFoam
    toolbox shared by Tian Tang, found at:
    https://bitbucket.org/tiantang/minigeotechfoam

    More reading at:

    T.Tang, O.Hededal & P.Cardiff (2014). On finite volume method implementation
    of poro-elasto-plasticity soil model. International Journal for
    Numerical and Analytical Methods in Geomechanics, 10.1002/nag.2361.

    T.Tang & O.Hededal (2014). Simulation of pore pressure accumulation under
    cyclic loading using finite volume method. Proceedings of the 8th European
    Conference on Numerical Methods in Geotechnical Engineering (NUMGE14),
    Volume 2, Pages 1301-1306.


Author
    Tian Tang, DTU.  All rights reserved.
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    poroLinGeomSolid.C

\*---------------------------------------------------------------------------*/

#ifndef poroLinGeomSolid_H
#define poroLinGeomSolid_H

#include "linSolid.H"
#include "totalDispSolid.H"
#include "SolverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class poroLinGeomSolid Declaration
\*---------------------------------------------------------------------------*/

class poroLinGeomSolid
:
    public linSolid<totalDispSolid>
{
    // Private data

        //- Poro pressure field
        volScalarField p_;

        //- Cell pressure gradient
        volVectorField gradp_;

        //- Hydraulic conductivity
        //  Note: Tang et al (2014) refers to this term as permeability
        dimensionedScalar hydraulicConductivity_;

        //- Specific weight of water
        //  Note: Tang et al (2014) refers to this term as water density
        dimensionedScalar gammaWater_;

        //- Porosity
        dimensionedScalar porosity_;

        //- Degree of saturation
        dimensionedScalar saturation_;

        //- Bulk modulus of water
        dimensionedScalar KWater_;

        //- Reciprocal of the effective bulk modulus of water and air mixture
        dimensionedScalar rKprime_;

    // Private Member Functions

        //- Check if the momentum and pressure equations are converged
        bool converged
        (
            const int iCorr,
            const SolverPerformance<vector>& solverPerfD,
            const SolverPerformance<scalar>& solverPerfp,
            const volVectorField& D,
            const volScalarField& p
        );

        //- Disallow default bitwise copy construct
        poroLinGeomSolid(const poroLinGeomSolid&);

        //- Disallow default bitwise assignment
        void operator=(const poroLinGeomSolid&);


public:

    //- Runtime type information
    TypeName("poroLinear");

    // Constructors

        //- Construct from components
        poroLinGeomSolid(dynamicFvMesh& mesh);


    // Destructor

        virtual ~poroLinGeomSolid()
        {}


    // Member Functions

        // Edit

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

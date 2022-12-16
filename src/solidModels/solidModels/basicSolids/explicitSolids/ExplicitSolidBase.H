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
    ExplicitSolidBase

Description
    Base class for explicit solid solvers

SourceFiles
    ExplicitSolidBase.C

\*---------------------------------------------------------------------------*/

#ifndef ExplicitSolidBase_H
#define ExplicitSolidBase_H

#include "mechanicalEnergies.H"
#include "dynamicRelaxation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                    Class ExplicitSolidBase Declaration
\*---------------------------------------------------------------------------*/

template<class IncrementalSolid>
class ExplicitSolidBase
:
    public IncrementalSolid
{
protected:
    // Protected data

        //- Wave speed field (compressive)
        surfaceScalarField wavespeed_;

        //- Wave speed field (shear)
        surfaceScalarField sWavespeed_;

        //- Energy calculation
        mechanicalEnergies energies_;

        //- Acceleration
        volVectorField a_;

        //- Dynamic relaxation
        dynamicRelaxation relaxation_;


    // Protected Member Functions

        //- Relax velocity field
        virtual void relax();

        //- Correct velocity boundary conditions
        virtual void correctUBCs(volVectorField& U);

        //- Return the divergence of stress
        virtual tmp<volVectorField> divStress() const = 0;

        //- Return the divergence of stress
        virtual tmp<surfaceVectorField> tractionSf() const = 0;

        //- Return first Piola tensor on a patch
        virtual tmp<vectorField> Pn(const fvPatch& patch) const = 0;


public:

    // Constructors

        //- Construct from components
        ExplicitSolidBase
        (
            const word& type,
            dynamicFvMesh& mesh,
            const nonLinearGeometry::nonLinearType nonLinear,
            const bool isSolid
        );

    // Destructor

        virtual ~ExplicitSolidBase()
        {}


    // Member Functions

        // Access

            //- Return the compressive wavespeed
            virtual const surfaceScalarField& wavespeed() const
            {
                return wavespeed_;
            }

            //- Return the shear wavespeed
            virtual const surfaceScalarField& sWavespeed() const
            {
                return sWavespeed_;
            }


        // Edit

            //- Solve the momentum equation
            void solveMomentum();

            //- Return the current Courant number
            virtual scalar CoNum() const;

            //- Return the max Courant number
            virtual scalar maxCoNum() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ExplicitSolidBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

Class
    Foam::twoPhaseCompressibleSystem

Description
    Compressible system contataining two equation of states.

SourceFiles
    twoPhaseCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseCompressibleSystem_H
#define twoPhaseCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseCompressibleSystem.H"
#include "twoPhaseFluidBlastThermo.H"
#include "fluxScheme.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class twoPhaseCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseCompressibleSystem
:
    public phaseCompressibleSystem
{
// Private data

    //- Primitative variables

        //- Equation of state
        twoPhaseFluidBlastThermo thermo_;

        //- Volume fraction of phase 1
        volScalarField& volumeFraction_;

        //- Phase 1 density
        volScalarField& rho1_;

        //- Phase 2 density
        volScalarField& rho2_;


    //- Conserved variables

        //- Phase 1 mass
        volScalarField alphaRho1_;

        //- phase 2 mass
        volScalarField alphaRho2_;


    //- Fluxes

        //- Volume fraction flux
        surfaceScalarField alphaPhi_;

        //- Phase 1 mass flux
        surfaceScalarField alphaRhoPhi1_;

        //- Phase 2 mass flux
        surfaceScalarField alphaRhoPhi2_;


        //- Calculate new alpha and rho fields
        virtual void calcAlphaAndRho();


public:

    TypeName("twoPhaseCompressibleSystem");

    // Constructor
    twoPhaseCompressibleSystem(const fvMesh& mesh);


    //- Destructor
    ~twoPhaseCompressibleSystem();


    // Member functions

        //- Decode primative variables
        void decode();

        //- Encode conserved variables
        void encode();

        //- Update fluxes
        void update();

        //- Return energy source term
        virtual tmp<volScalarField> ESource() const;

        //- Solve sub-step stepi
        virtual void solve();


    // Member Access Functions

        //- Const access to thermo model
        virtual const fluidBlastThermo& thermo() const
        {
            return thermo_;
        }

        //- Non-const access to thermo model
        virtual fluidBlastThermo& thermo()
        {
            return thermo_;
        }

        //- Return speed of sound
        virtual tmp<volScalarField> speedOfSound() const
        {
            return thermo_.speedOfSound();
        }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

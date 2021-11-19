/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::interfaceTwoPhaseCompressibleSystem

Description
    Compressible system contataining two equation of states that uses and
    interface.

SourceFiles
    interfaceTwoPhaseCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceTwoPhaseCompressibleSystem_H
#define interfaceTwoPhaseCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "compressibleBlastSystem.H"
#include "twoPhaseFluidBlastThermo.H"
#include "interfaceProperties.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceTwoPhaseCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class interfaceTwoPhaseCompressibleSystem
:
    public compressibleBlastSystem
{
// Private data

    //- Primitative variables

        //- Equation of state
        twoPhaseFluidBlastThermo& thermo_;

        //- Volume fraction of phase 1
        volScalarField& alpha1_;

        //- Phase 1 density
        volScalarField& rho1_;

        //- Volume fraction of phase 2
        volScalarField& alpha2_;

        //- Phase 2 density
        volScalarField& rho2_;

        //- Interface thickness parameter
        dimensionedScalar epsilon_;

        //- Level set function
        volScalarField levelSet_;


    //- Conserved variables

        //- Phase 1 mass
        volScalarField alphaRho1_;

        //- phase 2 mass
        volScalarField alphaRho2_;


    //- Fluxes

        //- Volume fraction flux for phase 1
        surfaceScalarField alphaPhi1_;

        //- Volume fraction flux for phase 2
        surfaceScalarField alphaPhi2_;

        //- Phase 1 mass flux
        surfaceScalarField alphaRhoPhi1_;

        //- Phase 2 mass flux
        surfaceScalarField alphaRhoPhi2_;

    //- Surface tension model
    interfaceProperties interfaceProperties_;

    //- Smooth distance function
    tmp<volScalarField> Psi(const volScalarField& alpha) const;

    //- Gradient of smooth distance function
    tmp<volVectorField> GradPhi(const volScalarField& alpha) const;

    //- Distance function
    tmp<volScalarField> D(const volScalarField& alpha) const;

    //- Interface compression function
    tmp<volScalarField> H() const;

    //- Realizable volume fraction sharpening
    tmp<volScalarField> sharpen(const volScalarField& alpha) const;


public:

    TypeName("interface");

    // Constructor
    interfaceTwoPhaseCompressibleSystem(const fvMesh& mesh);


    //- Destructor
    ~interfaceTwoPhaseCompressibleSystem();


    // Member functions

        //- Decode primative variables
        void decode();

        //- Encode conserved variables
        void encode();

        //- Update fluxes
        void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Solve implicit terms
        virtual void postUpdate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

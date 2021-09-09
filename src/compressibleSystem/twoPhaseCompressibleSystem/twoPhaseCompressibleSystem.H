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
    Foam::twoPhaseCompressibleSystem

Description
    Compressible system contataining two equation of states.

SourceFiles
    twoPhaseCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseCompressibleSystem_H
#define twoPhaseCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "compressibleBlastSystem.H"
#include "twoPhaseFluidBlastThermo.H"
#include "interfaceProperties.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class twoPhaseCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseCompressibleSystem
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

    //- Surface tension model
    autoPtr<interfaceProperties> interfacePropertiesPtr_;

    //- Allocate models
    virtual void setModels();


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

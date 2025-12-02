/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
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
    Foam::multiphaseCompressibleSystem

Description
    Compressible system containing more than two equation of states.

SourceFiles
    multiphaseCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseCompressibleSystem_H
#define multiphaseCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "compressibleBlastSystem.H"
#include "multiphaseFluidBlastThermo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class multiphaseCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class multiphaseCompressibleSystem
:
    public compressibleBlastSystem
{
protected:
// Protected data

    //- Primitative variables

        //- Equation of states
        multiphaseFluidBlastThermo& thermo_;

        //- Volume fractions
        PtrList<volScalarField>& alphas_;

        //- Phase densities
        UPtrList<volScalarField>& rhos_;


    //- Conserved variables

        //- Phase masses
        PtrList<volScalarField> alphaRhos_;


    //- Fluxes

        //- Volume fraction fluxes
        PtrList<surfaceScalarField> alphaPhis_;

        //- Phase mass fluxes
        PtrList<surfaceScalarField> alphaRhoPhis_;


    //- Cached changes from advection
    List<tmp<volScalarField>> alphaAdvection_;
    List<tmp<volScalarField>> alphaRhoAdvection_;


    //- Transport phase densities
    bool transportPhaseDensity_;

    //- Calculate reconstructed phaseDensity (alphaRho) with density
    bool densityReconstruction_;

    //- Use MULES volume fraction flux limiting
    bool MUSLESLimiting_;

    //- Phase to not modify volume fraction of
    labelList adjustablePhases_;

    //- Sharpen volume fraction fields
    List<bool> sharpen_;

        //- Calculate new alpha and rho fields
        virtual void calcAlphas();

        // Update fluxes given interpolated alphas and alphaRhos
        void updateFluxes
        (
            const PtrList<surfaceScalarField>& alphasOwn,
            const PtrList<surfaceScalarField>& alphasNei,

            const PtrList<surfaceScalarField>& alphaRhosOwn,
            const PtrList<surfaceScalarField>& alphaRhosNei
        );

        //- Return the index of the phase with the given name
        label whichPhase(const word& phaseName) const;

        //- Solve for density and volume fraction fields
        virtual void solveMass();


public:

    TypeName("multiphase");

    // Constructor
    multiphaseCompressibleSystem
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const bool initialize = true
    );


    //- Destructor
    ~multiphaseCompressibleSystem();


    // Member functions

        //- Decode primative variables
        virtual void decode();

        //- Encode conserved variables
        virtual void encode();

        //- Update fluxes
        virtual void update();

        //- Solve implicit terms
        virtual void postUpdate();

        //- Store changes from advection
        virtual void storeFluxDeltas();

        //- Clear temporary variables
        virtual void clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

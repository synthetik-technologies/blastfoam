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


        //- Calculate new alpha and rho fields
        virtual void calcAlphas();


public:

    TypeName("multiphaseCompressibleSystem");

    // Constructor
    multiphaseCompressibleSystem(const fvMesh& mesh);

    // Constructor
    multiphaseCompressibleSystem(const fvMesh& mesh, const bool);


    //- Destructor
    ~multiphaseCompressibleSystem();


    // Member functions

        //- Decode primative variables
        virtual void decode();

        //- Encode conserved variables
        virtual void encode();

        //- Update fluxes
        virtual void update();

        //- Solve a sub-step
        virtual void solve();

        //- Solve implicit terms
        virtual void postUpdate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

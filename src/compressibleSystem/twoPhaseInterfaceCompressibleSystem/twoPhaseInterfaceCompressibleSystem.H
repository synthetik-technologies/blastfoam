/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2022
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
    Foam::twoPhaseInterfaceCompressibleSystem

Description
    Compressible system contataining two equation of states.

SourceFiles
    twoPhaseInterfaceCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseInterfaceCompressibleSystem_H
#define twoPhaseInterfaceCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoPhaseCompressibleSystem.H"
#include "interfaceSystem.H"
#include "interfaceProperties.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class twoPhaseInterfaceCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseInterfaceCompressibleSystem
:
    public twoPhaseCompressibleSystem,
    public interfaceSystem
{
// Private data

    //- Interface properties
    interfaceProperties interface_;

    //- Level set function
    volScalarField psi_;

    //- Interface normal
    surfaceScalarField nHatf_;

    //- Default dictionary for GFM solve
    static dictionary GFMSolverDict_;
    //
    // //- Surface tension
    // autoPtr<surfaceTensionModel> surfaceTensionModel_;


    //- Momentum source
    virtual tmp<volVectorField> rhoUSource() const;

    //- Total energy source
    virtual tmp<volScalarField> rhoESource() const;


public:

    TypeName("twoPhaseInterface");

    // Constructor
    twoPhaseInterfaceCompressibleSystem(const fvMesh& mesh);


    //- Destructor
    ~twoPhaseInterfaceCompressibleSystem();


    // Member functions

        //- Decode primative variables
        void decode();

        //- Solve
        void solve();

        //- Update fluxes
        void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

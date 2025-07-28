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
#include "surfaceTensionModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class twoPhaseInterfaceCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseInterfaceCompressibleSystem
:
    public twoPhaseCompressibleSystem
{
// Private data

    //- Surface tension
    autoPtr<surfaceTensionModel> surfaceTensionModel_;

    //- Stabilisation for normalisation of the interface normal
    dimensionedScalar deltaN_;

    //- Interface normal
    surfaceScalarField nHatf_;


    // Protected member functions

        //- Add source terms to momentum and energy equations
        virtual void addSources
        (
            volVectorField::Internal& rhoUSource,
            volScalarField::Internal& rhoESource
        ) const;


public:

    TypeName("twoPhaseInterface");

    // Constructor
    twoPhaseInterfaceCompressibleSystem
    (
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    ~twoPhaseInterfaceCompressibleSystem();


    // Member functions

        // //- Decode primative variables
        // void decode();

        //- Update fluxes
        void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

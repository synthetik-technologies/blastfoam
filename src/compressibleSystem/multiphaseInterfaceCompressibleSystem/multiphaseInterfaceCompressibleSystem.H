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
    Foam::multiphaseInterfaceCompressibleSystem

Description
    Compressible system containing more than two equation of states.

SourceFiles
    multiphaseInterfaceCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseInterfaceCompressibleSystem_H
#define multiphaseInterfaceCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiphaseCompressibleSystem.H"
#include "multiphaseFluidBlastThermo.H"
#include "interfaceSystem.H"
#include "surfaceTensionModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class multiphaseInterfaceCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class multiphaseInterfaceCompressibleSystem
:
    public multiphaseCompressibleSystem,
    public interfaceSystem
{
    surfaceTensionTable stModels_;


public:

    TypeName("multiphaseInterface");

    // Constructor
    multiphaseInterfaceCompressibleSystem(const fvMesh& mesh);

    // Constructor
    multiphaseInterfaceCompressibleSystem(const fvMesh& mesh, const bool);


    //- Destructor
    ~multiphaseInterfaceCompressibleSystem();


    // Member functions

        //- Update fluxes
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

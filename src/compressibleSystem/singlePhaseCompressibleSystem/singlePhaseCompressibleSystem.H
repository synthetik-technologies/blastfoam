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
    Foam::singlePhaseCompressibleSystem

Description
    Compressible system with one phase

SourceFiles
    singlePhaseCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseCompressibleSystem_H
#define singlePhaseCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "compressibleBlastSystem.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class singlePhaseCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseCompressibleSystem
:
    public compressibleBlastSystem
{
public:

    TypeName("singlePhaseCompressibleSystem");

    // Constructor
    singlePhaseCompressibleSystem(const fvMesh& mesh);

    //- Destructor
    virtual ~singlePhaseCompressibleSystem();


    // Member functions

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

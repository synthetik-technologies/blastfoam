/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technology
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::reconstruction

Description
    Simple class for picking name of schemes and fields


SourceFiles
    ReconstructionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef reconstruction_H
#define reconstruction_H

#include "autoPtr.H"
#include "typeInfo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                 Class reconstruction Declaration
\*---------------------------------------------------------------------------*/

class reconstruction
{
public:

    //- Runtime type information
    TypeName("reconstruction");

    //- Constructor
    reconstruction();

    //- Destructor
    virtual ~reconstruction();


    // Static Functions

        //- Return the reconstruction scheme name
        static word scheme(const word&);
        static word scheme
        (
            const word&,
            const fvMesh&,
            const bool fail,
            const bool overwrite
        );
        static word scheme
        (
            const word&,
            const word&,
            const fvMesh&,
            const bool fail,
            const bool overwrite
        );

        //- Return the name of the owner field
        static word ownName(const word&);

        //- Return the name of the neighbour field
        static word neiName(const word&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

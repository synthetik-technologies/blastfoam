/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Renamed class
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
    Foam::blastSpecie

Description

SourceFiles
    blastSpecie.C

\*---------------------------------------------------------------------------*/

#ifndef specieBlast_H
#define specieBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"

using namespace Foam::constant::thermodynamic;

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class blastSpecie Declaration
\*---------------------------------------------------------------------------*/

class blastSpecie
{
protected:

    // Protected data

        //- Molecular weight
        scalar W_;


public:

    //- Return the instantiated type name
        static word typeName()
        {
            return word("specie");
        }

    // Constructor
    blastSpecie(const dictionary& dict)
    :
        W_(dict.subDict("specie").lookupType<scalar>("molWeight"))
    {}


    // Member Functions

        //- Return molecular weight
        scalar W() const
        {
            return W_;
        }

        //- Return molecular weight
        scalar R() const
        {
            return RR/W_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "blastSpecie.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

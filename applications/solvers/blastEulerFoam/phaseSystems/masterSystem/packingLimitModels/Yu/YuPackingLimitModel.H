/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::kineticTheoryModels::packingLimitModels::Yu

Description
TODO
    Reference:
    \verbatim
        Yu, A. B. and Standish, N. (1987).
        "Porosity calculations of multi-component mixtures of spherical
        particles."
        Powder Technology 52(3), 233-241
    \endverbatim

SourceFiles
    YuPackingLimitModel.C

\*---------------------------------------------------------------------------*/

#ifndef YuPackingLimitModel_H
#define YuPackingLimitModel_H

#include "binaryPackingLimitModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace packingLimitModels
{

/*---------------------------------------------------------------------------*\
                        Class Yu Declaration
\*---------------------------------------------------------------------------*/

class Yu
:
    public binary
{
protected:
// Protected member functions

    //- Interaction functions
    virtual scalar aij(const scalar rij) const;
    virtual scalar bij(const scalar rij) const;

public:

    //- Runtime type information
    TypeName("Yu");


    // Constructors

        //- Construct from components
        Yu
        (
            const dictionary& dict,
            const masterSystem& system
        );


    //- Destructor
    virtual ~Yu();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace packingLimitModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

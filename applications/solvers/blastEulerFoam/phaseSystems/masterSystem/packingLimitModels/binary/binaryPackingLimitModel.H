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
    Foam::kineticTheoryModels::packingLimitModels::YuStandish

Description

    Reference:
    \verbatim
        Yu, A. B. and Standish, N. (1987).
        "Porosity calculations of multi-component mixtures of spherical
        particles."
        Powder Technology 52(3), 233-241
    \endverbatim

SourceFiles
    binaryPackingLimitModel.C

\*---------------------------------------------------------------------------*/

#ifndef binaryPackingLimitModel_H
#define binaryPackingLimitModel_H

#include "packingLimitModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace packingLimitModels
{

/*---------------------------------------------------------------------------*\
                        Class YuStandish Declaration
\*---------------------------------------------------------------------------*/

class binary
:
    public packingLimitModel
{
protected:
// Protected member functions

    //- Interaction functions
    virtual scalar aij(const scalar rij) const = 0;
    virtual scalar bij(const scalar rij) const = 0;

public:

    //- Runtime type information
    TypeName("binary");


    // Constructors

        //- Construct from components
        binary
        (
            const dictionary& dict,
            const masterSystem& system
        );


    //- Destructor
    virtual ~binary();


    // Member functions

        //- Return packing limit
        virtual scalar alphaMax
        (
            const label celli,
            SortableList<scalar>& ds,
            const bool fixedDs
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace packingLimitModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

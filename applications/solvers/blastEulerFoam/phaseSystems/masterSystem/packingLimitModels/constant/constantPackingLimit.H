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
    Foam::packingLimitModels::constant

Description
    Constant packing limit

SourceFiles
    constantPackingLimit.C

\*---------------------------------------------------------------------------*/

#ifndef constantPackingLimit_H
#define constantPackingLimit_H

#include "packingLimitModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace packingLimitModels
{

/*---------------------------------------------------------------------------*\
               Class constantPackingLimit Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public packingLimitModel
{
    // Private data

        //- Maximum particle volume fraction
        mutable scalar maxAlpha_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from components
        constant
        (
            const dictionary& dict,
            const masterSystem& system
        );


    //- Destructor
    virtual ~constant();


    // Member functions

        //- Update the maximum packing limit
        virtual void updateAlphaMax
        (
            volScalarField::Internal& alphaMaxI
        ) const;

        //- Return packing limit
        virtual scalar alphaMax
        (
            const label celli,
            const SortableList<scalar>& ds
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace packingLimitModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

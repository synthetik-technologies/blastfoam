/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::pairInelastic

Description
    Pair forces between particles colliding with a spring, slider, damper model

\*---------------------------------------------------------------------------*/

#ifndef immersedpairInelastic_H
#define immersedpairInelastic_H

#include "pairModel.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class pairInelastic Declaration
\*---------------------------------------------------------------------------*/

class pairInelastic
:
    public pairModel
{
    // Private Data

public:

    //- Runtime type information
    TypeName("pairInelastic");


    // Constructors

        //- Construct from dictionary
        pairInelastic(const dictionary& dict);


    //- Destructor
    virtual ~pairInelastic();


    // Member Functions

        //- Calculate the pair interaction between parcels
        virtual void evaluatePair
        (
            immersedBoundaryObject& OA,
            immersedBoundaryObject& OB,
            const pairCollisionData& map
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

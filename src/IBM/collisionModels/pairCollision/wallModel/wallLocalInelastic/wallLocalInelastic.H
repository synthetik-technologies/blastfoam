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
    Foam::wallLocalInelastic

Description
    Inelastic collision between walls

\*---------------------------------------------------------------------------*/

#ifndef immersedwallLocalInelastic_H
#define immersedwallLocalInelastic_H

#include "wallModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class wallLocalInelastic Declaration
\*---------------------------------------------------------------------------*/

class wallLocalInelastic
:
    public wallModel
{
    // Private Data

        //- Coefficient of restitution
        scalar e_;




public:

    //- Runtime type information
    TypeName("wallLocalInelastic");


    // Constructors

        //- Construct from dictionary
        wallLocalInelastic(const dictionary& dict);


    //- Destructor
    virtual ~wallLocalInelastic();


    // Member Functions

        //- Calculate the wall interaction for a parcel
        virtual void evaluateWall
        (
            const polyPatch& patch,
            immersedBoundaryObject& object,
            const wallCollisionData& map
        ) const;

        //- Calculate the wall interaction for a parcel
        virtual void evaluateExternalForce
        (
            const polyPatch& patch,
            immersedBoundaryObject& object,
            const wallCollisionData& map
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

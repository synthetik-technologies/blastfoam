/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
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
    Foam::diameterReactionRates::length

Description
    Wrapper class to hold the surface reaction model for a diameter reaction
    rate

SourceFiles
    lengthDiameterReactionRate.C

\*---------------------------------------------------------------------------*/

#ifndef lengthDiameterReactionRate_H
#define lengthDiameterReactionRate_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "diameterReactionRate.H"
#include "surfaceReactionRate.H"
#include "runTimeSelectionTables.H"

namespace Foam
{
namespace diameterReactionRates
{

/*---------------------------------------------------------------------------*\
                    Class diameterReactionRate Declaration
\*---------------------------------------------------------------------------*/

class length
:
    public diameterReactionRate
{
protected:
// Protected data

    //- Reference to the diameter model
    autoPtr<surfaceReactionRate> rate_;

public:

    //- Runtime type information
    TypeName("length");


    // Constructor
    length
    (
        const diameterModel& dModel,
        autoPtr<surfaceReactionRate> rate
    );


    //- Destructor
    virtual ~length();


    // Public functions

        //- Reaction rate
        virtual tmp<volScalarField> dDdt
        (
            const volScalarField& p,
            const volScalarField& T
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterReactionRates
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

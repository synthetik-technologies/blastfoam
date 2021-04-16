/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
    Foam::reactionRates::tabulatedPRateActivation

Description
    tabulatedP reaction rate

SourceFiles
    tabulatedPReactionRate.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedPReactionRate_H
#define tabulatedPReactionRate_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reactionRate.H"
#include "lookupTable1D.H"

namespace Foam
{
namespace reactionRates
{

/*---------------------------------------------------------------------------*\
                    Class tabulatedP Declaration
\*---------------------------------------------------------------------------*/

class tabulatedP
:
    public reactionRate
{
// Private data

    //- lookup table
    lookupTable1D table_;

    //- Pressure scaling coefficient
    scalar pScale_;


public:

    //- Runtime type information
    TypeName("tabulatedP");

    // Constructor
    tabulatedP(const dictionary& dict);


    //- Destructor
    virtual ~tabulatedP();

    //- Construct and return a clone
    virtual autoPtr<reactionRate> clone() const
    {
        const tabulatedP& tprr(*this);
        return autoPtr<reactionRate>(new tabulatedP(tprr));
    }


    // Public functions

        //- Cell based reaction rate
        virtual scalar k(const scalar& p, const scalar& T) const;

        //- Reaction rate
        virtual tmp<volScalarField> k
        (
            const volScalarField& p,
            const volScalarField& T
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionRates
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

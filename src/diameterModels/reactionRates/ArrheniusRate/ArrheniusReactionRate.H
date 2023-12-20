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
    Foam::reactionRates::ArrheniusRateActivation

Description
    Arrhenius reaction rate

SourceFiles
    ArrheniusReactionRate.C

\*---------------------------------------------------------------------------*/

#ifndef ArrheniusReactionRate_H
#define ArrheniusReactionRate_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reactionRate.H"

namespace Foam
{
namespace reactionRates
{

/*---------------------------------------------------------------------------*\
                    Class Arrhenius Declaration
\*---------------------------------------------------------------------------*/

class Arrhenius
:
    public reactionRate
{
// Private data

    //- Temperature exponent
    dimensionedScalar beta_;
    
    //- Rate constant
    dimensionedScalar A_;

    //- Activation temperature
    dimensionedScalar Ta_;


public:

    //- Runtime type information
    TypeName("Arrhenius");

    // Constructor
    Arrhenius(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~Arrhenius();

    //- Construct and return a clone
    virtual autoPtr<reactionRate> clone() const
    {
        const Arrhenius& arr(*this);
        return autoPtr<reactionRate>(new Arrhenius(arr));
    }


    // Public functions

        //- Cell based reaction rate
        virtual scalar k
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const;

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

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
    Foam::surfaceReactionRates::ArrheniusRateActivation

Description
    Arrhenius reaction rate

SourceFiles
    ArrheniusSurfaceReactionRate.C

\*---------------------------------------------------------------------------*/

#ifndef ArrheniusSurfaceReactionRate_H
#define ArrheniusSurfaceReactionRate_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceReactionRate.H"

namespace Foam
{
namespace surfaceReactionRates
{

/*---------------------------------------------------------------------------*\
                    Class Arrhenius Declaration
\*---------------------------------------------------------------------------*/

class Arrhenius
:
    public surfaceReactionRate
{
// Private data

    //- Rate constant
    dimensionedScalar A_;

    //- Temperature exponent
    dimensionedScalar beta_;

    //- Activation temperature
    dimensionedScalar Ta_;


public:

    //- Runtime type information
    TypeName("Arrhenius");


    // Constructors

        //- Construct from dictionary
        Arrhenius(const dictionary& dict);

        //- Construct from runTime and dictionary
        Arrhenius(const Time& runTime, const dictionary& dict);

        //- Construct from fvMesh and dictionary
        Arrhenius(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~Arrhenius();

    //- Construct and return a clone
    virtual autoPtr<surfaceReactionRate> clone() const
    {
        const Arrhenius& arr(*this);
        return autoPtr<surfaceReactionRate>(new Arrhenius(arr));
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

        //- Reaction rate for a patch
        virtual tmp<scalarField> k
        (
            const fvPatchScalarField& p,
            const fvPatchScalarField& T
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceReactionRates
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

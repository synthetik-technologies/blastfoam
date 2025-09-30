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
    Foam::diameterReactionRates::HSR

Description
    Mass transfer rate based on arrhenius rate reaction for

    \verbatim
        Gurevich, M.A., Ozerova, G.E., Stepanov, A.M. (1970).
        Heterogeneous ignition of an aluminum particle in oxygen and water
        vapor.
        Combust Explos Shock Waves, 6 (3), 29-297.
    \endverbatim

SourceFiles
    HSRDiameterReactionRate.C

\*---------------------------------------------------------------------------*/

#ifndef HSRDiameterReactionRate_H
#define HSRDiameterReactionRate_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "diameterReactionRate.H"

namespace Foam
{
namespace diameterReactionRates
{

/*---------------------------------------------------------------------------*\
                    Class HSR Declaration
\*---------------------------------------------------------------------------*/

class HSR
:
    public diameterReactionRate
{
// Private data

    //- Rate constant
    dimensionedScalar A_;

    //- Activation temperature
    dimensionedScalar Ta_;

    //- Name ox oxidation specie
    const word oxidantName_;


public:

    //- Runtime type information
    TypeName("HSR");


    // Constructors

        //- Construct from runTime and dictionary
        HSR(const diameterModel& dModel, const dictionary& dict);


    //- Destructor
    virtual ~HSR();



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

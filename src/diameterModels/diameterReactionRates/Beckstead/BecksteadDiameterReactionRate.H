/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::diameterReactionRates::Beckstead

Description
    Burning rate of aluminium particles from the burn-time correlation of
    Beckstead (2005),

        t_b [ms] = a d0^n/(X_eff p^0.1 T^0.2),
        X_eff = X_O2 + 0.6 X_H2O + 0.22 X_CO2,

    with d0 in um, p in atm, T the gas temperature in K and X the mole
    fractions of the oxidisers in the gas. It is applied as the d^n law

        d(d^n)/dt = -X_eff p^0.1 T^0.2/a,

    i.e. dd/dt = -X_eff p^0.1 T^0.2/(a n d^(n-1)), so that a particle burns
    out in t_b at constant conditions. A particle ignites when its
    temperature reaches Tign; once its diameter has decreased below d0 it
    keeps burning (the diameter is transported with the particles, so this
    flag follows them), also when it is cooled below Tign.

    The oxidiser mole fractions are computed from the mass fractions of the
    oxidiser phase and the molecular weights of all its species, read from
    its physicalProperties dictionary; oxidisers that are not species of
    the phase are ignored.

    Reference:
    \verbatim
        Beckstead, M.W. (2005). Correlating aluminum burning times.
        Combustion, Explosion and Shock Waves 41(5), 533-546.
    \endverbatim

    Usage (in reactingDiameterCoeffs)
    \verbatim
        reactionRate    Beckstead;
        oxidiserPhase   gas;
        Tign            2000;       // Ignition temperature [K]
        a               0.00735;    // Optional, default 0.00735
        n               1.8;        // Optional, default 1.8
        oxidisers                   // Optional, default as below
        {
            O2      1;
            H2O     0.6;
            CO2     0.22;
        }
    \endverbatim

SourceFiles
    BecksteadDiameterReactionRate.C

\*---------------------------------------------------------------------------*/

#ifndef BecksteadDiameterReactionRate_H
#define BecksteadDiameterReactionRate_H

#include "diameterReactionRate.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterReactionRates
{

/*---------------------------------------------------------------------------*\
                          Class Beckstead Declaration
\*---------------------------------------------------------------------------*/

class Beckstead
:
    public diameterReactionRate
{
    // Private data

        //- Burn-time coefficient [ms um^-n]
        const scalar a_;

        //- Diameter exponent
        const scalar n_;

        //- Ignition temperature [K]
        const scalar Tign_;

        //- Name of the phase carrying the oxidisers
        const word oxidiserPhase_;

        //- Oxidiser species and their weights in X_eff
        HashTable<scalar> oxidisers_;

        //- Names of all species of the oxidiser phase
        wordList species_;

        //- Molecular weights of all species of the oxidiser phase
        scalarList W_;


    // Private Member Functions

        //- Read the species and their molecular weights
        void readSpecies();


public:

    //- Runtime type information
    TypeName("Beckstead");


    // Constructors

        //- Construct from the diameter model and a dictionary
        Beckstead
        (
            const diameterModel& dModel,
            const dictionary& dict
        );


    //- Destructor
    virtual ~Beckstead();


    // Member Functions

        //- Effective oxidiser mole fraction
        tmp<volScalarField> Xeff() const;

        //- Rate of change of the particle radius (positive) [m/s]
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

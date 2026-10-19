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
    flag follows them), also when it is cooled below Tign. With
    ignitedFraction (reactingDiameterModel) the ignition is tracked by the
    diameter model instead, and this is the rate of burning particles.

    The oxidiser mole fractions are computed from the mass fractions of the
    oxidiser phase and the molecular weights of all its species, read from
    its physicalProperties dictionary; oxidisers that are not species of
    the phase are ignored.

    Products (optional, products true): the burnt particle material is
    converted with the oxidisers of the cell to the products of a global
    reaction per oxidiser (default 4Al + 3O2 -> 2Al2O3, 2Al + 3H2O -> Al2O3
    + 3H2, 2Al + 3CO2 -> Al2O3 + 3CO), in the proportions of their terms in
    X_eff. productSpecies() and productYield() give the gas species and
    their mass per unit mass burnt (oxidisers negative); a reactingParticle
    mass transfer with productsFromReactionRate uses them, and the heat of
    reaction follows from the formation enthalpies of the species. The
    oxide is the solid species below the oxide melting temperature of the
    gas and the liquid species above it (if both are species of the gas).
    The rate is limited so that no more than oxidiserLimit of an oxidiser
    in the cell is consumed in one time step. Burning particles extinguish
    where the gas is colder than Textinction (default 0: never).

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
        products        true;       // Optional, default false
        oxidiserLimit   0.5;        // Optional, default 0.5
        Textinction     0;          // Optional, default 0
        oxideSolid      Al2O3(s);   // Optional, default Al2O3(s)
        oxideLiquid     Al2O3(L);   // Optional, default Al2O3(L)
        TmOxide         2327;       // Optional, default 2327
    \endverbatim
    An oxidiser entry can also give its weight and reaction: moles of
    particle material (Al) and of the gas species (consumed negative), e.g.
    \verbatim
            H2O
            {
                weight      0.6;
                Al          2;
                reaction    ((H2O -3) (Al2O3(s) 1) (H2 3));
            }
    \endverbatim

SourceFiles
    BecksteadDiameterReactionRate.C

\*---------------------------------------------------------------------------*/

#ifndef BecksteadDiameterReactionRate_H
#define BecksteadDiameterReactionRate_H

#include "diameterReactionRate.H"
#include "HashTable.H"
#include "Tuple2.H"

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

        //- Burn to the reaction products
        const bool products_;

        //- Largest fraction of an oxidiser consumed in one time step
        const scalar oxidiserLimit_;

        //- Gas temperature below which burning particles extinguish
        const scalar Textinction_;

        //- Ignition tracked by the diameter model (ignitedFraction)
        const bool externalIgnition_;

        //- Solid and liquid oxide species and the oxide melting temperature
        const word oxideSolid_;
        const word oxideLiquid_;
        const scalar TmOxide_;

        //- Reactions as read: per oxidiser the moles of particle material
        //  and of the gas species
        HashTable<scalar> nuP_;
        HashTable<List<Tuple2<word, scalar>>> reactions_;

        //- Mass of each gas species per unit mass burnt, per oxidiser
        HashTable<HashTable<scalar>> yields_;

        //- Gas species produced or consumed
        wordList productSpecies_;

        //- Is the liquid oxide a species of the gas
        bool liquidOxide_;


    // Private Member Functions

        //- Read the species and their molecular weights, and set the yields
        void readSpecies();

        //- Shares of the oxidisers in X_eff (sum to one where X_eff > 0)
        void shares
        (
            PtrList<volScalarField>& w,
            wordList& names
        ) const;


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

        //- Gas species produced or consumed (empty unless products)
        virtual wordList productSpecies() const;

        //- Mass of a species per unit mass burnt
        virtual tmp<volScalarField> productYield(const word& name) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterReactionRates
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

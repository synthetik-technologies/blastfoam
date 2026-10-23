/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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
    Foam::massTransferModels::reactingParticleMassTransfer

Description
    The change in volume of the diameterModel is used to determine mass
    transfer. Particles are assumed to be incompresssible.

    The species of the transferred mass are given as fractions
    (productSpecies, or <phase>Species), or, with
    productsFromReactionRate true, taken from the reaction rate of the
    reacting diameter model (e.g. Beckstead with products true): the
    products of the burning and the consumed oxidisers (negative), in
    yields that vary with the local oxidiser composition.

SourceFiles
    reactingParticleMassTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef reactingParticleMassTransfer_H
#define reactingParticleMassTransfer_H

#include "massTransferModel.H"
#include "phasePair.H"
#include "reactingDiameterModel.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace massTransferModels
{

/*---------------------------------------------------------------------------*\
                   Class reactingParticleMassTransfer Declaration
\*---------------------------------------------------------------------------*/

class reactingParticleMassTransfer
:
    public massTransferModel
{
// Private data

    //- Reacting phase name
    word reactingPhase_;

    //- Flag if the first phase is the reacting phase
    bool reactingIs1_;

    //- Name of species to transfer mass from
    mutable hashedWordList phase1Species_;

    //- Fraction of the phase1 species produced
    HashTable<scalar, word> phase1Yi_;

    //- Name of phase2 species
    mutable hashedWordList phase2Species_;

    //- Fraction of the product species produced
    HashTable<scalar, word> phase2Yi_;

    //- Take the product species and yields from the reaction rate
    bool productsFromRate_;

    //- Have the product species been taken from the reaction rate
    mutable bool productsSet_;

    //- Return the reacting phase model
    inline const phaseModel& reactingPhase() const
    {
        return reactingIs1_ ? pair_.phase1() : pair_.phase2();
    }

    //- Take the product species from the reaction rate (the diameter
    //  model is constructed after the mass transfer models)
    void setProductSpecies() const;

public:

    //- Runtime type information
    TypeName("reactingParticle");


    // Constructors

        //- Construct from components
        reactingParticleMassTransfer
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~reactingParticleMassTransfer();


    // Member Functions

        //- Return the mass transfer rate
        virtual tmp<volScalarField> K() const;

        //- Return the phase1 species of a given phase
        virtual tmp<volScalarField> phase1Y(const word& name) const;

        //- Return the phase1 species of a given phase
        virtual tmp<volScalarField> phase2Y(const word& name) const;

        //- Return the names of the dispersed phase mass fractions
        virtual hashedWordList phase1Species() const;

        virtual hashedWordList phase2Species() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace massTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::massTransferModels::reactingParticleMassTransfer

Description
    One phase is used to determine the interfacial Velocity

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

    //- Name of species to transfer mass from
    const hashedWordList reactingSpecies_;

    //- Fraction of the reacting species produced
    HashTable<scalar, word> reactingYi_;

    //- Name of product species
    const hashedWordList productSpecies_;

    //- Fraction of the product species produced
    HashTable<scalar, word> productYi_;


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

        //- Return the dispersed phase mass fraction of a given phase
        virtual tmp<volScalarField> dispersedYi(const word& name) const;

        //- Return the continuous phase mass fraction of a given phase
        virtual tmp<volScalarField> continuousYi(const word& name) const;

        //- Return the names of the dispersed phase mass fractions
        virtual hashedWordList dispersedSpecies() const
        {
            return reactingSpecies_;
        }

        //- Return the names of the continuous phase mass fractions
        virtual hashedWordList continuousSpecies() const
        {
            return productSpecies_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace massTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::massTransferModels::cavitation

Description
    The change in volume of the diameterModel is used to determine mass
    transfer. Particles are assumed to be incompresssible.

SourceFiles
    cavitationMassTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef cavitationMassTransfer_H
#define cavitationMassTransfer_H

#include "massTransferModel.H"
#include "phasePair.H"
#include "cavitationModel.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace massTransferModels
{

/*---------------------------------------------------------------------------*\
                   Class cavitation Declaration
\*---------------------------------------------------------------------------*/

class cavitation
:
    public massTransferModel
{
    //- Cavitation model
    autoPtr<cavitationModel> cavitation_;

public:

    //- Runtime type information
    TypeName("cavitation");


    // Constructors

        //- Construct from components
        cavitation
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~cavitation();


    // Member Functions

        //- Return the mass transfer rate
        virtual tmp<volScalarField> K() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace massTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

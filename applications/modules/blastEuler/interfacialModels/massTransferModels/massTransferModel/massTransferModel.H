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
    Foam::massTransferModel

Description

SourceFiles
    massTransferModel.C
    massTransferModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef massTransferModel_H
#define massTransferModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseModel.H"
#include "runTimeSelectionTables.H"
#include "hashedWordList.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                        Class massTransferModel Declaration
\*---------------------------------------------------------------------------*/

class massTransferModel
:
    public regIOobject
{
protected:

    // Protected data

        //- Reference to list of phaseModels
        const phasePair& pair_;


public:

    //- Runtime type information
    TypeName("massTransferModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            massTransferModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );


    // Static data members


    // Constructors

        //- Construct froma dictionary and a phase pair
        massTransferModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~massTransferModel();


    // Selectors

        static autoPtr<massTransferModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        const phasePair& pair() const
        {
            return pair_;
        }

        //- Return the mass transfer
        virtual tmp<volScalarField> K() const = 0;

        //- Return the species of a given phase
        tmp<volScalarField> Y(const phaseModel& phase, const word& name) const;

        //- Return the phase1 species of a given phase
        virtual tmp<volScalarField> phase1Y(const word& name) const;

        //- Return the phase1 species of a given phase
        virtual tmp<volScalarField> phase2Y(const word& name) const;

        //- Return the names of the given species
        virtual hashedWordList species(const phaseModel& phase) const;

        //- Return the names of the dispersed phase species
        virtual hashedWordList phase1Species() const
        {
            return hashedWordList::null();
        }

        //- Return the names of the continuous phase species
        virtual hashedWordList phase2Species() const
        {
            return hashedWordList::null();
        }

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const
        {
            return os.good();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

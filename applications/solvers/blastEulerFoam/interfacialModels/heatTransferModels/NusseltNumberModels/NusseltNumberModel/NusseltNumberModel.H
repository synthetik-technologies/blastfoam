/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technologies
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::NusseltNumberModel

Description

SourceFiles
    NusseltNumberModel.C

\*---------------------------------------------------------------------------*/

#ifndef NusseltNumberModel_H
#define NusseltNumberModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                        Class NusseltNumberModel Declaration
\*---------------------------------------------------------------------------*/

class NusseltNumberModel
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;


public:

    //- Runtime type information
    TypeName("NusseltNumberModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            NusseltNumberModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );


    // Static data members

        //- Coefficient dimensions
        static const dimensionSet dimK;


    // Constructors

        //- Construct from a dictionary and a phase pair
        NusseltNumberModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~NusseltNumberModel();


    // Selectors

        static autoPtr<NusseltNumberModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Nusselt number
        virtual tmp<volScalarField> Nu
        (
            const label nodei = -1,
            const label nodej = -1
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

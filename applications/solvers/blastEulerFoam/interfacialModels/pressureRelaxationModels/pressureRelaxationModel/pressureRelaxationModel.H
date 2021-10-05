/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
2017-05-24 Jeff Heylmun:    Added return functions for acceleration
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
    Foam::pressureRelaxationModel

Description

SourceFiles
    pressureRelaxationModel.C
    newpressureRelaxationModel.C

\*---------------------------------------------------------------------------*/

#ifndef pressureRelaxationModel_H
#define pressureRelaxationModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;
class swarmCorrection;

/*---------------------------------------------------------------------------*\
                           Class pressureRelaxationModel Declaration
\*---------------------------------------------------------------------------*/

class pressureRelaxationModel
:
    public regIOobject
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;


public:

    //- Runtime type information
    TypeName("pressureRelaxationModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            pressureRelaxationModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair,
                const bool registerObject
            ),
            (dict, pair, registerObject)
        );


    // Static data members

        //- Force coefficient dimensions
        static const dimensionSet dimMu;


    // Constructors

        // Construct without residual constants
        pressureRelaxationModel
        (
            const phasePair& pair,
            const bool registerObject
        );

        // Construct with residual constants
        pressureRelaxationModel
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~pressureRelaxationModel();


    // Selectors

        static autoPtr<pressureRelaxationModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Return the pressure relaxation coefficient
        virtual tmp<volScalarField> K
        (
            const label nodei = -1,
            const label nodej = -1
        ) const = 0;

        //- Return the pressure relaxation coefficient for a cell
        virtual scalar cellK
        (
            const label celli,
            const label nodei = -1,
            const label nodej = -1
        ) const = 0;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

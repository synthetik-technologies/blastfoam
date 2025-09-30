/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Syntyhetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivitive work of OpenFOAM.

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
    Foam::pressureRelaxationModels::constant

Description
    Constant rate pressure relaxation

SourceFiles
    constantPressureRelaxationModel.C

\*---------------------------------------------------------------------------*/

#ifndef constantPressureRelaxationModel_H
#define constantPressureRelaxationModel_H

#include "pressureRelaxationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace pressureRelaxationModels
{

/*---------------------------------------------------------------------------*\
                           Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public pressureRelaxationModel
{
    // Private data

        //- Pressure relaxation coefficient
        const dimensionedScalar mu_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from components
        constant
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~constant();


    // Member Functions

        //- Return the pressure relaxation coefficient
        virtual tmp<volScalarField> K
        (
            const label nodei = -1,
            const label nodej = -1
        ) const;

        //- Return the pressure relaxation coefficient for a cell
        virtual scalar cellK
        (
            const label celli,
            const label nodei = -1,
            const label nodej = -1
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pressureRelaxationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

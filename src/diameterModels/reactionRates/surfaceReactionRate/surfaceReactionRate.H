/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
    Foam::surfaceReactionRate

Description

SourceFiles
    surfaceReactionRate.C
    surfaceReactionRateNew.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceReactionRate_H
#define reactionRate_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class surfaceReactionRate Declaration
\*---------------------------------------------------------------------------*/

class surfaceReactionRate
{
public:

    //- Runtime type information
    TypeName("surfaceReactionRate");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceReactionRate,
            dictionary,
            (
                const Time& runTime,
                const dictionary& dict
            ),
            (runTime, dict)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceReactionRate,
            fvMesh,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );

    // Constructor
    surfaceReactionRate(const dictionary& dict);


    //- Destructor
    virtual ~surfaceReactionRate();


    // Selectors

        static autoPtr<surfaceReactionRate> New
        (
            const Time& runTime,
            const dictionary& dict
        );

        static autoPtr<surfaceReactionRate> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Construct and return a clone
    virtual autoPtr<surfaceReactionRate> clone() const = 0;


    // Public functions

        //- Optional update of necessary fields before calculating k
        virtual void update
        (
            const volScalarField& p,
            const volScalarField& T
        )
        {}

        //- Cell based reaction rate
        virtual scalar k
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const = 0;

        //- Reaction rate
        virtual tmp<volScalarField> k
        (
            const volScalarField& p,
            const volScalarField& T
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

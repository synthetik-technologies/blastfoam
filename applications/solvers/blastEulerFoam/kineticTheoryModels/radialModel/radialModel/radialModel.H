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
    Foam::kineticTheoryModels::radialModel

SourceFiles
    radialModel.C
    radialModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef radialModel_H
#define radialModel_H

#include "dictionary.H"
#include "volFields.H"
#include "kineticTheorySystem.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class radialModel Declaration
\*---------------------------------------------------------------------------*/

class radialModel
{
protected:

    // Protected data

        const dictionary& dict_;

        const kineticTheorySystem& kt_;


public:

    //- Runtime type information
    TypeName("radialModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        radialModel,
        dictionary,
        (
            const dictionary& dict,
            const kineticTheorySystem& kt
        ),
        (dict, kt)
    );


    // Constructors

        //- Construct from components
        radialModel
        (
            const dictionary& dict,
            const kineticTheorySystem& kt
        );

        //- Disallow default bitwise copy construction
        radialModel(const radialModel&) = delete;


    // Selectors

        static autoPtr<radialModel> New
        (
            const dictionary& dict,
            const kineticTheorySystem& kt
        );


    //- Destructor
    virtual ~radialModel();


    // Member Functions

        //- Radial distribution function
        virtual tmp<volScalarField> gs0
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const = 0;

        //- Derivative of the radial distribution function
        virtual tmp<volScalarField> gs0prime
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const = 0;

        virtual bool read()
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const radialModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

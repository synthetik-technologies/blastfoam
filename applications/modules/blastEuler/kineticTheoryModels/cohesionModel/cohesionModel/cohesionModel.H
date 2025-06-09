/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::kineticTheoryModels::cohesionModel

SourceFiles
    cohesionModel.C
    cohesionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef cohesionModel_H
#define cohesionModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class kineticTheoryModel;

namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                    Class cohesionModel Declaration
\*---------------------------------------------------------------------------*/

class cohesionModel
{
protected:

    // Protected data

        const kineticTheoryModel& kt_;


public:

    //- Runtime type information
    TypeName("cohesionModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        cohesionModel,
        dictionary,
        (
            const dictionary& dict,
            const kineticTheoryModel& kt
        ),
        (dict, kt)
    );


    // Constructors

        //- Construct from components
        cohesionModel
        (
            const dictionary& dict,
            const kineticTheoryModel& kt
        );

        //- Disallow default bitwise copy construction
        cohesionModel(const cohesionModel&) = delete;


    // Selectors

        static autoPtr<cohesionModel> New
        (
            const dictionary& dict,
            const kineticTheoryModel& kt
        );


    //- Destructor
    virtual ~cohesionModel();


    // Member Functions

        //- Viscosity
        virtual tmp<volScalarField> nu() const = 0;

        //- Granular pressure coefficient
        virtual tmp<volScalarField> Ps() const = 0;

        //- Derivative of the granular pressure w.r.t. volume  fraction
        virtual tmp<volScalarField> dPsdAlpha() const = 0;

        //- Derivative of the granular pressure w.r.t. granular temperature
        virtual tmp<volScalarField> dPsdTheta() const = 0;

        //- Dissipation source
        virtual tmp<volScalarField> dissipationSource
        (
            const dimensionedScalar&
        ) const = 0;


        //- Update saved fields
        virtual void update()
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cohesionModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

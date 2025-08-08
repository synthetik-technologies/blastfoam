/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::reactionRateFlameAreaModel

Description
    Abstract class for reaction rate per flame area unit

SourceFiles
    reactionRateFlameAreaModel.C
    reactionRateFlameAreaModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef reactionRateFlameAreaModel_H
#define reactionRateFlameAreaModel_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "volFields.H"
#include "blastCombustionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class fvMesh;

/*---------------------------------------------------------------------------*\
                    Class reactionRateFlameAreaModel Declaration
\*---------------------------------------------------------------------------*/

class reactionRateFlameAreaModel
{
protected:
    // Protected data

        //- Dictionary
        dictionary coeffDict_;

        //- Mesh reference
        const fvMesh& mesh_;

        //- Combustion model owner
        const blastCombustionModel& combModel_;

        //- Fuel name
        word fuel_;

        //- Fuel consumption rate per unit of flame area
        volScalarField omega_;

public:
            class consumptionSpeed
            {
                // Private Data

                    //- Maximum consumption speed
                    scalar omega0_;

                    //- Exponential factor
                    scalar eta_;

                    //- Extinction strain
                    scalar sigmaExt_;

                    //- Minimum consumption speed
                    scalar omegaMin_;


                // Private Member Functions

                    //- Return consumption rate
                    scalar omega0Sigma(scalar sigma, scalar a) const;


            public:

                // Constructors

                    //- Construct from dictionary
                    consumptionSpeed(const dictionary& dict);


                //- Destructor
                virtual ~consumptionSpeed();


                // Member Functions

                    //- Return speed consumption rate temp
                    tmp<volScalarField> omega0Sigma(const volScalarField& sigma);


                    // Access functions

                        inline scalar omega0() const
                        {
                            return omega0_;
                        }

                        inline scalar eta() const
                        {
                            return eta_;
                        }

                        inline scalar sigmaExt() const
                        {
                            return sigmaExt_;
                        }

                        inline scalar omegaMin() const
                        {
                            return omegaMin_;
                        }


                // IO

                    //- Update properties
                    void read(const dictionary& dict);
            };


public:

    //- Runtime type information
    TypeName("reactionRateFlameArea");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            reactionRateFlameAreaModel,
            dictionary,
            (
                const word modelType,
                const dictionary& dict,
                const fvMesh& mesh,
                const blastCombustionModel& combModel
            ),
            (modelType, dict, mesh, combModel)
        );


    // Constructors

        //- Construct from components
        reactionRateFlameAreaModel
        (
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh,
            const blastCombustionModel& combModel
        );

        //- Disallow default bitwise copy construction
        reactionRateFlameAreaModel(const reactionRateFlameAreaModel&) = delete;


    // Selector

        static autoPtr<reactionRateFlameAreaModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const blastCombustionModel& combModel
        );


    // Destructor

        virtual ~reactionRateFlameAreaModel();


    // Member Functions

        //- Access functions

            //- Return omega
            const volScalarField& omega() const
            {
                return omega_;
            }


        //- Correct omega
        virtual void correct(const volScalarField& sigma) = 0;

        //- Update from dictionary
        virtual bool read(const dictionary& dictProperties);


    // Member Operators

         //- Disallow default bitwise assignment
        void operator=(const reactionRateFlameAreaModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
20-06-2020 Synthetik Applied Technologies: | Time of arrival implementation
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
    Foam::functionObjects::blastQuantities

Description
    Calculates standard blast quantities including overpressure, impulse,
    time of arrival, and positive phase duration

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    blastQuantities.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_blastQuantities_H
#define functionObjects_blastQuantities_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class blastQuantities Declaration
\*---------------------------------------------------------------------------*/

class blastQuantities
:
    public fvMeshFunctionObject
{
public:

    // Enumeration of stages of Friedlander curve
    enum FriedlanderStage
    {
        PRE = -1,
        POSITIVE_PHASE = 0,
        NEGATIVE_PHASE = 1,
        POST = 2
    };

protected:
    // Private Member data

        //- Name of pressure field
        word pName_;

        //- Reference pressure
        dimensionedScalar pRef_;

        //- Overpressure field name
        word overpressureName_;

        //- Impulse field
        volScalarField impulse_;

        //- Time of arrival field
        volScalarField timeOfArrival_;

        //- Positive phase duratuion
        volScalarField positivePhaseDuration_;

        //- Negative phase duratuion
        volScalarField negativePhaseDuration_;

        //- Current stage of Friedlander fit
        volScalarField stage_;

        //- Max pressure field
        volScalarField& pMax_;

        //- Return the given field or create
        volScalarField& lookupOrCreate
        (
            const word& name,
            const dimensionSet& dims
        ) const;


        void update
        (
            const scalar p,
            scalar& pMax,
            scalar& stage,
            scalar& toa,
            scalar& ppd,
            scalar& npd
        )
        {
            const scalar op = p - pRef_.value();

            // Still havent hit time of arrival
            if (p > pMax)
            {
                toa = obr_.time().value();
                pMax = p;

                // Reset positive phase and negative phase durations
                ppd = 0.0;
                npd = 0.0;
            }
            // Positivive over pressure and in PRE or POSITIVE_PHASE
            else if
            (
                op > 0
            && stage <= POSITIVE_PHASE
            )
            {
                stage = POSITIVE_PHASE;
                ppd += obr_.time().deltaTValue();
            }
            // Negative overpressure, in POSITIVE_PHASE or NEGATIVE_PHASE
            else if
            (
                op < 0
             && (
                    stage == POSITIVE_PHASE
                 || stage == NEGATIVE_PHASE
                )
            )
            {
                stage = NEGATIVE_PHASE;
                npd += obr_.time().deltaTValue();
            }
            // Hit positive overpressure after negative phase, fitting to
            // friedlander curve is now invalid
            else if
            (
                op > 0
            && stage >= NEGATIVE_PHASE
            )
            {
                stage = POST;
            }
        }


public:

    //- Runtime type information
    TypeName("blastQuantities");

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        blastQuantities
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        blastQuantities(const blastQuantities&) = delete;


    //- Destructor
    virtual ~blastQuantities();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Calculate the field averages
        virtual bool execute();

        //- Write the field averages
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const blastQuantities&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

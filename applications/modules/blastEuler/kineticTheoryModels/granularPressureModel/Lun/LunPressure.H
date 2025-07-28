/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::Lun

Description

SourceFiles
    LunPressure.C

\*---------------------------------------------------------------------------*/

#ifndef LunPressure_H
#define LunPressure_H

#include "granularPressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace granularPressureModels
{

/*---------------------------------------------------------------------------*\
                           Class Lun Declaration
\*---------------------------------------------------------------------------*/

class Lun
:
    public granularPressureModel
{

public:

    //- Runtime type information
    TypeName("Lun");


    // Constructors

        //- Construct from components
        Lun(const dictionary& dict, const kineticTheorySystem& kt);


    //- Destructor
    virtual ~Lun();


    // Member Functions

        //- Granular pressure coefficient
        virtual tmp<volScalarField> granularPressure
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& Theta1,
            const volScalarField& Theta2,
            const volScalarField& g0,
            const dimensionedScalar& e
        ) const;

        //- Derivative of the granular pressure w.r.t. volume  fraction
        virtual tmp<volScalarField> granularPressureByAlpha
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& Theta1,
            const volScalarField& Theta2,
            const volScalarField& g0,
            const volScalarField& g0prime,
            const dimensionedScalar& e
        ) const;

        //- Derivative of the granular pressure w.r.t. granular temperature
        virtual tmp<volScalarField> granularPressureByTheta
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& Theta1,
            const volScalarField& Theta2,
            const volScalarField& g0,
            const dimensionedScalar& e
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace granularPressureModels
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::activationModels::ArrheniusRateActivation

Description
    Activation rate based on two step Arrhenius rate model

SourceFiles
    ArrheniusRateActivation.C

\*---------------------------------------------------------------------------*/

#ifndef ArrheniusRateActivation_H
#define ArrheniusRateActivation_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "activationModel.H"

namespace Foam
{
namespace activationModels
{

/*---------------------------------------------------------------------------*\
                           Class ArrheniusRateActivation Declaration
\*---------------------------------------------------------------------------*/

class ArrheniusRateActivation
:
    public activationModel
{
// Private data

    //- reference density
    dimensionedScalar rho0_;

    //- Diameter
    dimensionedScalar dp_;

    //- Switch temperature
    dimensionedScalar Ts_;


    // Low temperature coefficients

        //- Activation Rate
        dimensionedScalar ALow_;

        //- Activation energy
        dimensionedScalar TaLow_;


    // High temperature coefficients

        //- Activation Rate
        dimensionedScalar AHigh_;

        //- Activation energy
        dimensionedScalar TaHigh_;

        //- Reference to temperature
        const volScalarField& T_;


    //- Return the time rate of chage of lambda
    virtual tmp<volScalarField> delta() const;


public:

    //- Runtime type information
    TypeName("ArrheniusRate");

    // Constructor
    ArrheniusRateActivation
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );


    //- Destructor
    virtual ~ArrheniusRateActivation();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace activationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

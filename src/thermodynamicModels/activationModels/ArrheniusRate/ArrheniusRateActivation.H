/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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
#include "diameterModel.H"

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
    autoPtr<diameterModel> dModel_;

    //- Reference to diameter
    const volScalarField& dp_;

    //- Ignition temperature
    dimensionedScalar Tign_;

    //- Switch temperature
    dimensionedScalar Ts_;


    // Low temperature coefficients

        //- Activation Rate
        dimensionedScalar ALow_;

        //- Activation energy
        dimensionedScalar EaLow_;


    // High temperature coefficients

        //- Activation Rate
        dimensionedScalar AHigh_;

        //- Activation energy
        dimensionedScalar EaHigh_;

        //- Name of temperature field
        word TName_;


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


    // Public functions

        //- Update
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Post update
        virtual void postUpdate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace activationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
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
    Foam::turbulentDispersionModels::Panicker

Description
    Bubble dispersion model of Panicker et al. Used to describe the drag
    induced by one bubble on neighbour bubbles, and to achieve hyperbolicity
    of the two-phase model.

    References:
    \verbatim
        "On the hyperbolicity of two-fluid models for gas-liquid flows."
        Panicker, N., Passalacqua, A., Fox, R. O.,
        In preparation
    \endverbatim

SourceFiles
    Panicker.C

\*---------------------------------------------------------------------------*/

#ifndef Panicker_H
#define Panicker_H

#include "turbulentDispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace turbulentDispersionModels
{

/*---------------------------------------------------------------------------*\
                             Class Panicker Declaration
\*---------------------------------------------------------------------------*/

class Panicker
:
    public turbulentDispersionModel
{
    // Private data

        //- Schmidt number
        const dimensionedScalar Cdis_;

        //- Residual phase fraction
        const dimensionedScalar residualAlpha_;


public:

    //- Runtime type information
    TypeName("Panicker");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Panicker
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~Panicker();


    // Member Functions

        //- Turbulent diffusivity
        //  multiplying the gradient of the phase-fraction
        virtual tmp<volScalarField> D(const label, const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulentDispersionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

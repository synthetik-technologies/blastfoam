/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2025-06-09 Jeff Heylmun:    Added cell based returns
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
    Foam::heatTransferModels::constantCoefficient

Description
    Constant heat transfer coefficient, optionally scaled by a field (e.g.
    the reacted fraction lambda of an explosive phase, so that unreacted
    explosive does not equilibrate with the gas):

    \verbatim
        Kh      1e15;
        field   lambda.compBAl;     // optional
    \endverbatim

SourceFiles
    constantCoefficientHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef constantCoefficientHeatTransfer_H
#define constantCoefficientHeatTransfer_H

#include "dispersedHeatTransferModel.H"
#include "NuHeatTransfer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                          Class constantCoefficient Declaration
\*---------------------------------------------------------------------------*/

class constantCoefficient
:
    public dispersedHeatTransferModel
{

    // Heat transfer coefficient
    dimensionedScalar Kh_;

    //- Name of the field scaling the coefficient (optional)
    word fieldName_;

    //- Return the scaling field, or null if none
    const volScalarField* field() const;

public:

    //- Runtime type information
    TypeName("constantCoefficient");


    // Constructors

        //- Construct from components
        constantCoefficient
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~constantCoefficient();


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        virtual tmp<volScalarField> Ki() const;

        //- The heat transfer function K used in the enthalpy equation
        virtual scalar cellKi(const label celli) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

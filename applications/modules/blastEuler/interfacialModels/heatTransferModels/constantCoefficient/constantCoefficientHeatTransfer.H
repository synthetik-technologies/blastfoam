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
    Constant heat transfer coefficient

SourceFiles
    constantCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef constantCoefficient_H
#define constantCoefficient_H

#include "heatTransferModel.H"
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
    public heatTransferModel
{

    // Heat transfer coefficient
    dimensionedScalar Kh_;

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
        virtual tmp<volScalarField> K() const;

        //- The heat transfer function K used in the enthalpy equation
        virtual scalar cellK(const label celli) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

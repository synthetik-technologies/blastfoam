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
    Foam::heatTransferModels::noneNu

Description
    No heat transfer but calculates a Nusselt number

SourceFiles
    noneNu.C

\*---------------------------------------------------------------------------*/

#ifndef noneNu_H
#define noneNu_H

#include "heatTransferModel.H"
#include "NuHeatTransfer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                          Class noneNu Declaration
\*---------------------------------------------------------------------------*/

class noneNu
:
    public NuHeatTransfer
{

public:

    //- Runtime type information
    TypeName("noneNu");


    // Constructors

        //- Construct from components
        noneNu
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~noneNu();


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        tmp<volScalarField> K
        (
            const label nodei = -1,
            const label nodej = -1
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

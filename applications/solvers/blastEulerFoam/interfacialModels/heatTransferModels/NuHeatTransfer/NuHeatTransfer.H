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
    Foam::heatTransferModels::NuHeatTransfer

Description
    Base class for a heat transfer model that uses a Nusselt number

SourceFiles
    RanzMarshall.C

\*---------------------------------------------------------------------------*/

#ifndef NuHeatTransfer_H
#define NuHeatTransfer_H

#include "heatTransferModel.H"
#include "NusseltNumberModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                          Class NuHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class NuHeatTransfer
:
    public heatTransferModel
{
protected:
// Protected data

    //- Nusselt number model
    autoPtr<NusseltNumberModel> NuModel_;


public:

    //- Runtime type information
    TypeName("NuHeatTransfer");


    // Constructors

        //- Construct from components
        NuHeatTransfer
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~NuHeatTransfer();


    // Member Functions

        //- Function to return Nusselt number model
        const NusseltNumberModel& NusseltNo() const
        {
            return NuModel_();
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

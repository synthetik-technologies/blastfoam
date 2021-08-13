/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technologies
     \\/     M anipulation  |
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
    Foam::NusseltNumberModels::GelperinEpstein

Description
    Nusselt number correlation of GelperinEpstein

SourceFiles
    GelperinEpsteinNusseltNumberModel.C

\*---------------------------------------------------------------------------*/

#ifndef GelperinEpsteinNusseltNumberModel_H
#define GelperinEpsteinNusseltNumberModel_H

#include "NusseltNumberModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace NusseltNumberModels
{

/*---------------------------------------------------------------------------*\
                          Class GelperinEpstein Declaration
\*---------------------------------------------------------------------------*/

class GelperinEpstein
:
    public NusseltNumberModel
{
public:

    //- Runtime type information
    TypeName("GelperinEpstein");


    // Constructors

        //- Construct from components
        GelperinEpstein
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~GelperinEpstein();


    // Member Functions

        //- Nusselt number
        tmp<volScalarField> Nu
        (
            const label nodei = -1,
            const label nodej = -1
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace NusseltNumberModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

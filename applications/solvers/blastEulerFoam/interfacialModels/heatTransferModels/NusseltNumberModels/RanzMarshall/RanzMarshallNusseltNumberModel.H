/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022 Synthetik Applied Technologies
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
    Foam::NusseltNumberModels::RanzMarshall

Description
    Ranz Marshall Nusselt number for spheres

SourceFiles
    RanzMarshall.C

\*---------------------------------------------------------------------------*/

#ifndef RanzMarshallNusseltNumber_H
#define RanzMarshallNusseltNumber_H

#include "NusseltNumberModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace NusseltNumberModels
{

/*---------------------------------------------------------------------------*\
                          Class constant Declaration
\*---------------------------------------------------------------------------*/

class RanzMarshall
:
    public NusseltNumberModel
{
public:

    //- Runtime type information
    TypeName("RanzMarshall");


    // Constructors

        //- Construct from components
        RanzMarshall
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~RanzMarshall();


    // Member Functions

        //- Nusselt number
        virtual tmp<volScalarField> Nu
        (
            const label nodei = -1,
            const label nodej = -1
        ) const;

        //- Nusselt number for cell
        virtual scalar cellNu
        (
            const label,
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

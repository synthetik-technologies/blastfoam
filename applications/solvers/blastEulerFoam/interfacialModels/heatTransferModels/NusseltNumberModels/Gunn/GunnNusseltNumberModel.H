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
    Foam::NusseltNumberModels::Gunn

Description
    Nusselt number correlation of Gunn

SourceFiles
    GunnNusseltNumberModel.C

\*---------------------------------------------------------------------------*/

#ifndef GunnNusseltNumberModel_H
#define GunnNusseltNumberModel_H

#include "NusseltNumberModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace NusseltNumberModels
{

/*---------------------------------------------------------------------------*\
                          Class Gunn Declaration
\*---------------------------------------------------------------------------*/

class Gunn
:
    public NusseltNumberModel
{
public:

    //- Runtime type information
    TypeName("Gunn");


    // Constructors

        //- Construct from components
        Gunn
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~Gunn();


    // Member Functions

        //- Nusselt number
        tmp<volScalarField> Nu(const label nodei, const label nodej) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace NusseltNumberModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

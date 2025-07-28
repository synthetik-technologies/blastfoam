/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
2025-06-09 Jeff Heylmun:    Added cell based returns
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
    Foam::liftModels::constantLiftCoefficient

Description
    Constant coefficient lift model.

SourceFiles
    constantLiftCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef constantLiftCoefficient_H
#define constantLiftCoefficient_H

#include "dispersedLiftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{

/*---------------------------------------------------------------------------*\
                           Class constantLiftCoefficient Declaration
\*---------------------------------------------------------------------------*/

class constantLiftCoefficient
:
    public dispersedLiftModel
{
    // Private data

        //- Constant lift coefficient
        const dimensionedScalar Cl_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from a dictionary and a phase pair
        constantLiftCoefficient
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~constantLiftCoefficient();


    // Member Functions

        //- Lift coefficient
        virtual tmp<volScalarField> Cl() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace liftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

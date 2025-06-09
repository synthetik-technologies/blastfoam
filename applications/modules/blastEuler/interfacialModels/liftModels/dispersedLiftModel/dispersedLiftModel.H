/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
2017-05-24 Jeff Heylmun:    Added return functions for acceleration
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
    Foam::liftModels::dispersedLiftModel

Description

SourceFiles
    dispersedLiftModel.C

\*---------------------------------------------------------------------------*/

#ifndef dispersedLiftModel_H
#define dispersedLiftModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "liftModel.H"
#include "phasePair.H"

namespace Foam
{
namespace liftModels
{

/*---------------------------------------------------------------------------*\
                Class liftModels::dispersedLiftModel Declaration
\*---------------------------------------------------------------------------*/

class dispersedLiftModel
:
    public liftModel
{
public:

    //- Runtime type information
    TypeName("dispersedLiftModel");


    // Constructors

        // Construct with residual constants
        dispersedLiftModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~dispersedLiftModel();


    // Member Functions

        //- Lift coefficient
        virtual tmp<volScalarField> Cl() const = 0;

        //- Return the phase-intensive lift force
        virtual tmp<volVectorField> Fi() const;

        //- Return the lift force
        virtual tmp<volVectorField> F() const;

        //- Return the lift force
        virtual tmp<surfaceScalarField> Ff() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace liftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
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
    Foam::turbulentDispersionModels::dispersedTurbulentDispersionModel

Description

SourceFiles
    dispersedTurbulentDispersionModel.C

\*---------------------------------------------------------------------------*/

#ifndef dispersedTurbulentDispersionModel_H
#define dispersedTurbulentDispersionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "turbulentDispersionModel.H"
#include "phasePair.H"
#include "phaseCompressibleMomentumTransportModel.H"

namespace Foam
{
namespace turbulentDispersionModels
{

/*---------------------------------------------------------------------------*\
            Class dispersedTurbulentDispersionModel Declaration
\*---------------------------------------------------------------------------*/

class dispersedTurbulentDispersionModel
:
    public turbulentDispersionModel
{
public:

    //- Runtime type information
    TypeName("dispersedTurbulentDispersionModel");

    // Constructors

        //- Construct from a dictionary and a phase pair
        dispersedTurbulentDispersionModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~dispersedTurbulentDispersionModel();


    // Member Functions

        //- Return reference to continuous phase turbulence model
        const phaseCompressible::momentumTransportModel&
            continuousTurbulence() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulentDispersionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

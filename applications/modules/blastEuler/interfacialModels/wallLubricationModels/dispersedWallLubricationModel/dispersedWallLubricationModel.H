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
    Foam::dispersedWallLubricationModel

Description

SourceFiles
    dispersedWallLubricationModel.C
    dispersedWallLubricationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef dispersedWallLubricationModel_H
#define dispersedWallLubricationModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallLubricationModel.H"
#include "phasePair.H"

namespace Foam
{
namespace wallLubricationModels
{

/*---------------------------------------------------------------------------*\
                      Class dispersedWallLubricationModel Declaration
\*---------------------------------------------------------------------------*/

class dispersedWallLubricationModel
:
    public wallLubricationModel
{
public:

    //- Runtime type information
    TypeName("dispersedWallLubricationModel");

    // Constructors

        //- Construct from components
        dispersedWallLubricationModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~dispersedWallLubricationModel();


    // Member Functions

        //- Return phase-intensive wall lubrication force
        virtual tmp<volVectorField> Fi() const = 0;

        //- Return wall lubrication force
        virtual tmp<volVectorField> F() const;

        //- Return face wall lubrication force
        virtual tmp<surfaceScalarField> Ff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallLubricationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

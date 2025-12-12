/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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
    Foam::virtualMassModels::Zuber

Description
    Virtual mass model of Zuber.

    Reference:
    \verbatim
        "On the dispersed two-phase flow in the laminar regime"
        Zuber, N.,
        Chemical Engineering Science,
        issue 19, 897-917, 1964
    \endverbatim

SourceFiles
    ZuberVirtualMass.C

\*---------------------------------------------------------------------------*/

#ifndef ZuberVirtualMass_H
#define ZuberVirtualMass_H

#include "dispersedVirtualMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace virtualMassModels
{

/*---------------------------------------------------------------------------*\
                           Class Zuber Declaration
\*---------------------------------------------------------------------------*/

class Zuber
:
    public dispersedVirtualMassModel
{
// Private data

    //- Residual alpha
    dimensionedScalar residualAlpha_;


public:

    //- Runtime type information
    TypeName("Zuber");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Zuber
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~Zuber();


    // Member Functions

        //- Virtual mass coefficient
        virtual tmp<volScalarField> Cvm() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace virtualMassModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::dragModels::Stokes

Description
    Length based drag model.

SourceFiles
    lengthBased.C

\*---------------------------------------------------------------------------*/

#ifndef lengthBased_H
#define lengthBased_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class lengthBased Declaration
\*---------------------------------------------------------------------------*/

class lengthBased
:
    public dragModel
{
    // Private data

        //- C coefficient
        const dimensionedScalar C_;


public:

    //- Runtime type information
    TypeName("lengthBased");


    // Constructors

        //- Construct from components
        lengthBased
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~lengthBased();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe
        (
            const label nodei = -1,
            const label nodej = -1
        ) const;

        //- The drag function used in the momentum equation
        virtual tmp<volScalarField> KI
        (
            const label nodei = -1,
            const label nodej = -1
        ) const;

        //- Drag coefficient
        virtual scalar CdRei
        (
            const label,
            const label nodei = -1,
            const label nodej = -1
        ) const;

        //- Return the phase-intensive drag coefficient Ki
        //  used in the momentum equations
        //    ddt(alpha1*rho1*U1) + ... = ... alphad*K*(U1-U2)
        //    ddt(alpha2*rho2*U2) + ... = ... alphad*K*(U2-U1)
        virtual scalar KIi
        (
            const label,
            const label nodei = -1,
            const label nodej = -1
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

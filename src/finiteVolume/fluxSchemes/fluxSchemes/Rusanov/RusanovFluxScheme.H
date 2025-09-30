/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022 Synthetik Applied Technologies
     \\/     M anipulation  |
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
    Foam::fluxSchemes::Rusanov

Description

    References:
    \verbatim
        Rusanov, P. L.  (1981).
        Approximate Riemann solvers, parameter vectors, and difference schemes.
        Journal of Computational Physics, 43(2), 357-372.
    \endverbatim

SourceFiles
    RusanovFluxScheme.C

\*---------------------------------------------------------------------------*/

#ifndef RusanovFluxScheme_H
#define RusanovFluxScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluxScheme.H"

namespace Foam
{
namespace fluxSchemes
{

/*---------------------------------------------------------------------------*\
                           Class Rusanov Declaration
\*---------------------------------------------------------------------------*/

class Rusanov
:
    public fluxScheme
{

    // Saved variables
    tmp<surfaceScalarField> phivOwn_;
    tmp<surfaceScalarField> phivNei_;
    tmp<surfaceScalarField> lambda_;


    // Private functions

        //- Calcualte fluxes
        virtual void calculateFluxes
        (
            const scalar& rhoOwn, const scalar& rhoNei,
            const vector& UOwn, const vector& UNei,
            const scalar& eOwn, const scalar& eNei,
            const scalar& pOwn, const scalar& pNei,
            const scalar& cOwn, const scalar& cNei,
            const vector& Sf,
            scalar& phi,
            scalar& rhoPhi,
            vector& rhoUPhi,
            scalar& rhoEPhi,
            const label facei, const label patchi = -1
        );

        //- Calculate energy flux for an addition internal energy
        virtual scalar energyFlux
        (
            const scalar& rhoOwn, const scalar& rhoNei,
            const vector& UOwn, const vector& UNei,
            const scalar& eOwn, const scalar& eNei,
            const scalar& pOwn, const scalar& pNei,
            const vector& Sf,
            const label facei, const label patchi = -1
        ) const;

        //- Interpolate field
        virtual scalar interpolate
        (
            const scalar& fOwn, const scalar& fNei,
            const label facei, const label patchi = -1
        ) const;

        //- Calcualte flux given own/nei values
        virtual scalar calculateFlux
        (
            const scalar& fOwn, const scalar& fNei,
            const scalar& phi,
            const label facei, const label patchi = -1
        ) const;


public:

    //- Runtime type information
    TypeName("Rusanov");

    // Constructor
    Rusanov(const surfaceScalarField& phi);


    //- Destructor
    virtual ~Rusanov();


    // Member Functions

        //- Clear savedFields
        virtual void clear();

        //- Allocate saved fields
        virtual void createSavedFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluxSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

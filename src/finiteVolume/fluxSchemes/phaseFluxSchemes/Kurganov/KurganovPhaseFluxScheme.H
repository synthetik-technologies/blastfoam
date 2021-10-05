/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2019-10-21  Jeff Heylmun:   Moved from rhoCentralFoam to runtime selectable
                            method.
-------------------------------------------------------------------------------License
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
    Foam::phaseFluxSchemes::Kurganov

Description

SourceFiles
    KurganovPhaseFluxScheme.C

\*---------------------------------------------------------------------------*/

#ifndef KurganovPhaseFluxScheme_H
#define KurganovPhaseFluxScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseFluxScheme.H"

namespace Foam
{
namespace phaseFluxSchemes
{

/*---------------------------------------------------------------------------*\
                           Class Kurganov Declaration
\*---------------------------------------------------------------------------*/

class Kurganov
:
    public phaseFluxScheme
{

    // Saved variables

        tmp<surfaceScalarField> aPhivOwn_;
        tmp<surfaceScalarField> aPhivNei_;
        tmp<surfaceScalarField> aOwn_;
        tmp<surfaceScalarField> aNei_;
        tmp<surfaceScalarField> aSf_;


    // Private functions

        virtual void calculateFluxes
        (
            const scalar& alphaOwn, const scalar& alphaNei,
            const scalar& rhoOwn, const scalar& rhoNei,
            const vector& UOwn, const vector& UNei,
            const scalar& eOwn, const scalar& eNei,
            const scalar& pOwn, const scalar& pNei,
            const scalar& cOwn, const scalar& cNei,
            const vector& Sf,
            scalar& phi,
            scalar& alphaPhi,
            scalar& alphaRhoPhi,
            vector& alphaRhoUPhi,
            scalar& alphaRhoEPhi,
            const label facei, const label patchi = -1
        );

        //- Calcualte fluxes
        virtual void calculateFluxes
        (
            const scalar& alphaOwn, const scalar& alphaNei,
            const scalar& rhoOwn, const scalar& rhoNei,
            const scalarList& alphasOwn, const scalarList& alphasNei,
            const scalarList& rhosOwn, const scalarList& rhosNei,
            const vector& UOwn, const vector& UNei,
            const scalar& eOwn, const scalar& eNei,
            const scalar& pOwn, const scalar& pNei,
            const scalar& cOwn, const scalar& cNei,
            const vector& Sf,
            scalar& phi,
            scalarList& alphaPhis,
            scalarList& alphaRhoPhis,
            vector& alphaRhoUPhi,
            scalar& alphaRhoEPhi,
            const label facei, const label patchi = -1
        );

        //- Interpolate field
        virtual scalar interpolate
        (
            const scalar& fOwn, const scalar& fNei,
            const bool rho,
            const label facei, const label patchi = -1
        ) const;


public:

    //- Runtime type information
    TypeName("Kurganov");

    // Constructor
    Kurganov(const fvMesh& mesh, const word& name);


    //- Destructor
    virtual ~Kurganov();


    // Member Functions

        //- Clear savedFields
        virtual void clear();

        //- Allocate saved fields
        virtual void createSavedFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseFluxSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

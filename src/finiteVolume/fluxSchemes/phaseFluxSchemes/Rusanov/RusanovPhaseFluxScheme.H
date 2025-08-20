/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2019-10-21  Jeff Heylmun:   Moved from rhoCentralFoam to runtime selectable
                            method.
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
    Foam::phaseFluxSchemes::Rusanov

Description

SourceFiles
    RusanovPhaseFluxScheme.C

\*---------------------------------------------------------------------------*/

#ifndef RusanovPhaseFluxScheme_H
#define RusanovPhaseFluxScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseFluxScheme.H"

namespace Foam
{
namespace phaseFluxSchemes
{

/*---------------------------------------------------------------------------*\
                           Class Rusanov Declaration
\*---------------------------------------------------------------------------*/

class Rusanov
:
    public phaseFluxScheme
{
protected:
    // Saved variables

        tmp<surfaceScalarField> lambda_;


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
            scalar& alphaRhoPhi,
            vector& alphaRhoUPhi,
            scalar& alphaRhoEPhi,
            const label facei, const label patchi = -1
        );

        //- Calcualte interpolated alpha values
        //  U \cdot (\sum_i Sf_i alphaf_1) + \sum_i alphaf_corr,i
        virtual scalar calculateAlphaCorrector
        (
            const scalar& alphaOwn, const scalar& alphaNei,
            const label facei, const label patchi = -1
        ) const;

        //- Calcualte flux given own/nei values
        virtual scalar calculateFlux
        (
            const scalar& fOwn, const scalar& fNei,
            const scalar& phi,
            const label facei, const label patchi = -1
        ) const;

        //- Interpolate field
        virtual scalar interpolate
        (
            const scalar& fOwn, const scalar& fNei,
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

} // End namespace phaseFluxSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

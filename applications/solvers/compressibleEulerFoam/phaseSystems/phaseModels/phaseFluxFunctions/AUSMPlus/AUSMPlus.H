/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
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
    Foam::fluxSchemes::AUSMPlus

Description
    Computes flux using AUSM+ scheme

    References:
    \verbatim
        Liou, Meng-Sing (1996).
        A Sequel to AUSM: AUSM+.
        Journal of Computational Physics, 129, 364-382.
    \endverbatim

SourceFiles
    AUSMPlus.C

\*---------------------------------------------------------------------------*/

#ifndef AUSMPlus_H
#define AUSMPlus_H

#include "fluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxSchemes
{

/*---------------------------------------------------------------------------*\
                  Class AUSMPlus Declaration
\*---------------------------------------------------------------------------*/

class AUSMPlus
:
    public fluxScheme
{
    // Private Data

        //- Coefficients
        scalar ku_;
        scalar kp_;

        //- Cut off mach Number
        scalar cutOffMa_;

        //- Residual density
        scalar residualRho_;

        //- Minimum velocity
        scalar residualU_;

        //- Saved flux
        tmp<surfaceScalarField> phi_;


    // Private functions

        //- First polynomial
        scalar f
        (
            const scalar& M,
            const label sign
        ) const;

        //- Second polynomial
        scalar beta
        (
            const scalar& M,
            const label sign,
            const scalar& fa
        ) const;

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

        //- Interpolate field
        virtual scalar interpolate
        (
            const scalar& fOwn, const scalar& fNei,
            const bool rho,
            const label facei, const label patchi = -1
        ) const;


public:

    //- Runtime type information
    TypeName("AUSM+");

    // Costructor
    AUSMPlus(const fvMesh& mesh, const word& name);


    //- Destructor
    virtual ~AUSMPlus();


    // Member Functions

        //- Clear savedFields
        virtual void clear();

        //- Allocate saved fields
        virtual void createSavedFields();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace fluxSchemes

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

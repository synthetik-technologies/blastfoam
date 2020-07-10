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
    Foam::fluxSchemes::AUSMPlusUp

Description
    Computes flux using AUSM+up scheme with low speed correction

    References:
    \verbatim
        Liou, Meng-Sing (2006).
        A sequel to AUSM, Part II: AUSM+-up for all speeds.
        Journal of Computational Physics, 214, 137-170.
    \endverbatim

SourceFiles
    AUSMPlusUp.C

\*---------------------------------------------------------------------------*/

#ifndef AUSMPlusUp_H
#define AUSMPlusUp_H

#include "fluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxSchemes
{

/*---------------------------------------------------------------------------*\
                  Class AUSMPlusUp Declaration
\*---------------------------------------------------------------------------*/

class AUSMPlusUp
:
    public fluxScheme
{
    // Private Data

        //- Coefficients
        scalar beta_ = 0.125;
        scalar Kp_ = 0.25;
        scalar Ku_ = 0.75;
        scalar sigma_ = 1.0;
        scalar fa_ = 1.0;
        scalar alpha_ = 3.0/16.0*(5.0*sqr(fa_) - 4.0);

        //- Saved flux
        tmp<surfaceScalarField> phi_;


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

        //- Calcualte fluxes
        virtual void calculateFluxes
        (
            const scalarList& alphasOwn, const scalarList& alphasNei,
            const scalarList& rhosOwn, const scalarList& rhosNei,
            const scalar& rhoOwn, const scalar& rhoNei,
            const vector& UOwn, const vector& UNei,
            const scalar& eOwn, const scalar& eNei,
            const scalar& pOwn, const scalar& pNei,
            const scalar& cOwn, const scalar& cNei,
            const vector& Sf,
            scalar& phi,
            scalarList& alphaPhis,
            scalarList& alphaRhoPhis,
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
            const label facei, const label patchi = -1
        ) const;

        //- Interpolate field
        virtual scalar interpolate
        (
            const scalar& fOwn, const scalar& fNei,
            const label facei, const label patchi = -1
        ) const;

        // Mach number polynomials

            //- First order
            scalar M1(const scalar& Ma, const label sign) const
            {
                return 0.5*(Ma + sign*mag(Ma));
            }

            //- Second order
            scalar M2(const scalar& Ma, const label sign) const
            {
                if (mag(Ma) >= 1.0)
                {
                    return M1(Ma, sign);
                }
                return sign*0.25*sqr(Ma + sign*1.0);
            }

            //- Fourth order
            scalar M4(const scalar& Ma, const label sign) const
            {
                if (mag(Ma) >= 1.0)
                {
                    return M1(Ma, sign);
                }
                return sign*0.25*sqr(Ma + sign) + sign*beta_*sqr(sqr(Ma) - 1.0);
            }

            //- Fifth order
            scalar P5(const scalar& Ma, const label sign) const
            {
                if (mag(Ma) >= 1.0)
                {
                    return 0.5*(1.0 + sign*Foam::sign(Ma));
                }
                return
                    0.25*sqr(Ma + sign)*(2.0 - sign*Ma)
                  + sign*alpha_*Ma*sqr(sqr(Ma) - 1.0);
            }


public:

    //- Runtime type information
    TypeName("AUSM+up");

    // Costructor
    AUSMPlusUp(const fvMesh& mesh);


    //- Destructor
    virtual ~AUSMPlusUp();


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

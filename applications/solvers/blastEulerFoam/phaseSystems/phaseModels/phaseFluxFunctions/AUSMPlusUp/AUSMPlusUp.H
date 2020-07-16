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
    Computes flux using AUSM+Up scheme for a granular phase

    References:
    \verbatim
        Houim, R. W., Oran, E. S. (2017).
        A multiphase model for compressible granular-gaseous flows:
        formulation and initial tests,
        Journal of Fluid Mechanics, 789, 166-220
    \endverbatim

SourceFiles
    AUSMPlusUp.C

\*---------------------------------------------------------------------------*/

#ifndef AUSMPlusUp_H
#define AUSMPlusUp_H

#include "fluxScheme.H"
#include "kineticTheorySystem.H"

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
        scalar beta_;
        scalar fa_;
        scalar D_;

        //- Reference to packing limit
        const kineticTheorySystem& kineticTheorySystem_;

        //- Interpolated total volume fraction
        tmp<surfaceScalarField> alphapOwn_;
        tmp<surfaceScalarField> alphapNei_;

        //- Interpolated packing limit
        tmp<surfaceScalarField> alphaMaxf_;

        //- Interpolated min friction
        tmp<surfaceScalarField> alphaMinFrictionf_;

        //- Cut off mach Number
        scalar cutOffMa_;

        //- Residual phasic density
        scalar residualAlphaRho_;

        //- Minimum velocity
        scalar residualU_;

        //- Saved flux
        tmp<surfaceScalarField> phi_;


    // Private functions

        //- First polynomial
        scalar M1
        (
            const scalar& M,
            const label sign
        );

        //- Second polynomial
        scalar M2
        (
            const scalar& M,
            const label sign
        );

        //- 4th polynomial
        scalar M4
        (
            const scalar& M,
            const label sign,
            const scalar& beta
        );

        //- Fith polynomial
        scalar P5
        (
            const scalar& M,
            const label sign,
            const scalar& xi
        );

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
        )
        {
            NotImplemented;
        }

        //- Interpolate field
        virtual scalar interpolate
        (
            const scalar& fOwn, const scalar& fNei,
            const bool rho,
            const label facei, const label patchi = -1
        ) const;

        //- Update fields before calculating fluxes
        virtual void preUpdate(const volScalarField& p);

        //- Correct fluxes
        virtual void postUpdate();


public:

    //- Runtime type information
    TypeName("AUSM+Up");

    // Costructor
    AUSMPlusUp(const fvMesh& mesh, const word& name);


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

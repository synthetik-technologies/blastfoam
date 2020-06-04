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
    Foam::fluxSchemes::HLLCP

Description

    References:
    \verbatim
        Xie, W., Zheng, R., Lai, J., Li, H.  (2018).
        An accurate and robust HLLC-type Riemann solver for the compressible
        Euler system at various Mach numbers.
        Internation Journal of Numberical Methods for Fluids, 89, 430-463.
    \endverbatim

SourceFiles
    HLLCP.C

\*---------------------------------------------------------------------------*/

#ifndef HLLCP_H
#define HLLCP_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluxScheme.H"

namespace Foam
{
namespace fluxSchemes
{

/*---------------------------------------------------------------------------*\
                           Class HLLCP Declaration
\*---------------------------------------------------------------------------*/

class HLLCP
:
    public fluxScheme
{

    // Saved variables

        tmp<surfaceScalarField> SOwn_;
        tmp<surfaceScalarField> SNei_;
        tmp<surfaceScalarField> SStar_;
        tmp<surfaceScalarField> pStar_;
        tmp<surfaceScalarField> phip_;
        tmp<surfaceVectorField> UTilde_;
        tmp<surfaceScalarField> UvOwn_;
        tmp<surfaceScalarField> UvNei_;

        tmp<surfaceScalarField> f_;

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

        //- Update fields before calculating fluxes
        virtual void preUpdate(const volScalarField& p);

        //- Correct fluxes
        virtual void postUpdate();


public:

    //- Runtime type information
    TypeName("HLLCP");

    // Constructor
    HLLCP(const fvMesh& mesh);


    //- Destructor
    virtual ~HLLCP();


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

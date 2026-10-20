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
    Foam::phaseFluxSchemes::HLLC

Description
    HLLC flux of a phase. With hybridHLL (fluxSchemes.<phase> dictionary,
    default no) the mass, momentum and energy fluxes are blended towards
    the HLL flux near shocks, which removes the carbuncle of strong shocks
    aligned with the mesh. The HLL weight of a face is 1 - f^3, f the
    smallest pressure ratio (min(p_i/p_j, p_j/p_i)) between a cell of the
    face and its neighbours (the HLLCP sensor of Xie et al. 2018), with the
    mixture pressure where the trace treatment is on. Faces with a trace
    side keep the HLLC flux.

    References:
    \verbatim
        Toro, E.F., Spruce, M., Speares, W.  (1993).
        Restoration of the contact surface in the HLL-Riemann solver.
        Shock Waves, 4, 25-34.

        Xie, W., Zheng, R., Lai, J., Li, H.  (2018).
        An accurate and robust HLLC-type Riemann solver for the compressible
        Euler system at various Mach numbers.
        International Journal for Numerical Methods in Fluids, 89, 430-463.
    \endverbatim

SourceFiles
    HLLCPhaseFluxScheme.C

\*---------------------------------------------------------------------------*/

#ifndef HLLCPhaseFluxScheme_H
#define HLLCPhaseFluxScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseFluxScheme.H"

namespace Foam
{
namespace phaseFluxSchemes
{

/*---------------------------------------------------------------------------*\
                           Class HLLC Declaration
\*---------------------------------------------------------------------------*/

class HLLC
:
    public phaseFluxScheme
{

    // Saved variables

        tmp<surfaceScalarField> SOwn_;
        tmp<surfaceScalarField> SNei_;
        tmp<surfaceScalarField> SStar_;

        //- Blend towards HLL near shocks
        Switch hybridHLL_;

        //- HLL weight of the faces
        tmp<surfaceScalarField> thetaHLL_;

    // Private functions

        //- Set the HLL weight of the faces from the pressure
        virtual void preUpdate(const volScalarField& p);

        //- Clear the HLL weight
        virtual void postUpdate();

        //- Calcualte fluxes
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
        //  U \cdot (\sum_i Sf_i alphaf) + \sum_i alphaf_corr,i
        virtual scalar calculateAlphaCorrector
        (
            const scalar& alphaOwn, const scalar& alphaNei,
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
    TypeName("HLLC");

    // Constructor
    HLLC
    (
        const surfaceScalarField& phi,
        const scalar residualAlpha
    );


    //- Destructor
    virtual ~HLLC();


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

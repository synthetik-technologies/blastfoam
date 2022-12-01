/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::shallowWaterSystem

Description
    System of equations for shallow water equations

SourceFiles
    shallowWaterSystem.C

\*---------------------------------------------------------------------------*/

#ifndef shallowWaterSystem_H
#define shallowWaterSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegrationSystem.H"
#include "shallowFluxScheme.H"
#include "rainfallModel.H"
#include "uniformDimensionedFields.H"
#include "Function3Fwd.H"
#include "NamedEnum.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class shallowWaterSystem Declaration
\*---------------------------------------------------------------------------*/

class shallowWaterSystem
:
    public timeIntegrationSystem
{
public:
    // Public enumeration
    enum Friction
    {
        none,
        Manning,
        DarcyWeisbach
    };

    static const NamedEnum<Friction, 3> frictionTypes;

protected:
// Protected data

    //- Flood properties
    IOdictionary dict_;

    //- Gravitational acceleration
    uniformDimensionedVectorField g_;

    //- Gravitational acceleration magnitude
    dimensionedScalar magg_;

    //- Are rotational effects included
    bool rotating_;

    //- Angular velocity
    dimensionedVector omega_;

    //- Centripital acceleration
    dimensionedVector F_;


    //- Primitative variables

        //- Height
        volScalarField h_;

        //- Base height
        volScalarField h0_;

        //- Velocity
        volVectorField U_;

        //- Height-momentum
        volVectorField hU_;


    //- Fluxes

        //- Volumetric flux
        surfaceScalarField phi_;

        //- Height flux
        surfaceScalarField hPhi_;

        //- Hight momentum flux
        surfaceVectorField hUPhi_;


    //- Is friction included
    bool friction_;
    Friction frictionType_;

    //- Manning's coefficient (friction)
    autoPtr<volScalarField> nPtr_;

    //- Coefficient of friction
    autoPtr<volScalarField> fPtr_;

    //- Viscous dissipation
    bool viscous_;

    //- Turbulent diffusion
    bool turbulence_;

    //- Horizontal eddy viscosity
    dimensionedScalar muh_;

    //- Vertical eddy viscosity
    dimensionedScalar muv_;

    //- Laminar friction coefficient
    dimensionedScalar kl_;

    //- Turbulent friction coefficient
    dimensionedScalar kt_;

    //- Is rain included
    bool rain_;

    //- Rainfall model
    autoPtr<rainfallModel> rainfall_;

    //- Momentum source (geometric)
    autoPtr<Function3<vector>> S0Ptr_;


    //- Flux scheme
    autoPtr<shallowFluxScheme> flux_;

    //- Minimum height
    dimensionedScalar hMin_;


public:

    TypeName("shallow");

    // Constructor
    shallowWaterSystem(const fvMesh& mesh, const bool read = true);

    //- Destructor
    virtual ~shallowWaterSystem();


    // Member functions

        //- Decode primative variables
        virtual void decode();

        //- Encode conserved variables
        virtual void encode();

        //- Update fluxes
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Post update
        virtual void postUpdate();


    // Member Access Functions

        //- Return height
        const volScalarField& h() const
        {
            return h_;
        }

        //- Return base height
        const volScalarField& h0() const
        {
            return h0_;
        }

        //- Return velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return volumetric flux
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Return the Courant number
        scalar CoNum() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

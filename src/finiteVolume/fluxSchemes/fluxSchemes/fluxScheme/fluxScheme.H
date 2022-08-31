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
    Foam::fluxScheme

Description
    Base class for flux schemes to interpolate fields and loop over faces
    and boundaries for a single shared velocity and energy

SourceFiles
    fluxScheme.C
    fluxSchemeNew.C

\*---------------------------------------------------------------------------*/

#ifndef fluxScheme_H
#define fluxScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluxSchemeBase.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluxScheme Declaration
\*---------------------------------------------------------------------------*/

class fluxScheme
:
    public fluxSchemeBase
{
protected:

    //- Dictionary for fluxScheme
    const dictionary& dict_;

    //- Saved interpolated U field
    tmp<surfaceVectorField> Uf_;


    // Protected Functions

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
        ) = 0;

        //- Update
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
            scalarList& alphaRhosPhis,
            vector& rhoUPhi,
            scalar& rhoEPhi,
            const label facei, const label patchi = -1
        ) = 0;

        //- Calculate energy flux for an addition internal energy
        virtual scalar energyFlux
        (
            const scalar& rhoOwn, const scalar& rhoNei,
            const vector& UOwn, const vector& UNei,
            const scalar& eOwn, const scalar& eNei,
            const scalar& pOwn, const scalar& pNei,
            const vector& Sf,
            const label facei, const label patchi = -1
        ) const = 0;

        //- Update fields before calculating fluxes
        virtual void preUpdate(const volScalarField& p)
        {}

        //- Correct fluxes
        virtual void postUpdate()
        {}


public:

    //- Runtime type information
    TypeName("fluxScheme");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluxScheme,
            singlePhase,
            (const fvMesh& mesh),
            (mesh)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluxScheme,
            multiphase,
            (const fvMesh& mesh),
            (mesh)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluxScheme,
            interface,
            (const fvMesh& mesh),
            (mesh)
        );

    // Constructor
    fluxScheme(const fvMesh& mesh);


    //- Destructor
    virtual ~fluxScheme();


    // Selectors

        static autoPtr<fluxScheme> NewSingle(const fvMesh& mesh);
        static autoPtr<fluxScheme> NewMulti(const fvMesh& mesh);
        static autoPtr<fluxScheme> NewInterface(const fvMesh& mesh);


    // Member Functions

        //- Clear savedFields
        virtual void clear();

        //- Allocate saved fields
        virtual void createSavedFields();

        //- Single phase flux
        void update
        (
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& e,
            const volScalarField& p,
            const volScalarField& c,
            surfaceScalarField& phi,
            surfaceScalarField& rhoPhi,
            surfaceVectorField& rhoUPhi,
            surfaceScalarField& rhoEPhi
        );

        //- Multiphase fluxes
        void update
        (
            const PtrList<volScalarField>& alphas,
            const UPtrList<volScalarField>& rhos,
            const volVectorField& U,
            const volScalarField& e,
            const volScalarField& p,
            const volScalarField& c,
            surfaceScalarField& phi,
            PtrList<surfaceScalarField>& alphaPhis,
            PtrList<surfaceScalarField>& alphaRhoPhis,
            surfaceScalarField& rhoPhi,
            surfaceVectorField& rhoUPhi,
            surfaceScalarField& rhoEPhi
        );

        //- Two phase fluxes
        void update
        (
            const volScalarField& alpha,
            const volScalarField& rho1,
            const volScalarField& rho2,
            const volVectorField& U,
            const volScalarField& e,
            const volScalarField& p,
            const volScalarField& c,
            surfaceScalarField& phi,
            surfaceScalarField& alphaPhi,
            surfaceScalarField& alphaRhoPhi1,
            surfaceScalarField& alphaRhoPhi2,
            surfaceScalarField& alphaRhoPhi,
            surfaceVectorField& rhoUPhi,
            surfaceScalarField& rhoEPhi
        );

        //- Calculate energy flux for an addition internal energy
        tmp<surfaceScalarField> energyFlux
        (
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& e,
            const volScalarField& p
        ) const;

        //- Return interpolated U field
        tmp<surfaceVectorField> Uf() const;

        //- Return Anti-diffustion coefficient for a given volume fraction
        //  and update anti-diffusion flux
        virtual tmp<surfaceScalarField> snGradAlpha
        (
            const volScalarField& alpha
        ) const;

        //- Return Anti-diffustion coefficient for a given volume fraction
        virtual tmp<surfaceScalarField> AD
        (
            const volScalarField& alpha
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
public:

    //- Return the flux scheme
    static const fluxScheme& findFluxScheme
    (
        const surfaceScalarField& phi
    )
    {
        return dynamicCast<const fluxScheme>
        (
            fluxSchemeBase::findFluxScheme(phi)
        );
    }

protected:

    //- Dictionary for fluxScheme
    const dictionary& dict_;

    //- Saved interpolated U field
    mutable tmp<surfaceVectorField> Uf_;


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


        //- Return the name of a give field
        word fieldName(const word& name) const
        {
            return type() + ":" + name;
        }


public:

    //- Runtime type information
    TypeName("fluxScheme");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluxScheme,
            singlePhase,
            (const surfaceScalarField& phi),
            (phi)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluxScheme,
            multiphase,
            (const surfaceScalarField& phi),
            (phi)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluxScheme,
            interface,
            (const surfaceScalarField& phi),
            (phi)
        );

    // Constructor
    fluxScheme(const surfaceScalarField& phi);


    //- Destructor
    virtual ~fluxScheme();


    // Selectors

        static autoPtr<fluxScheme> NewSingle(const surfaceScalarField& phi);
        static autoPtr<fluxScheme> NewMulti(const surfaceScalarField& phi);
        static autoPtr<fluxScheme> NewInterface(const surfaceScalarField& phi);


    // Member Functions

        //- Clear savedFields
        virtual void clear();

        //- Allocate saved fields
        virtual void createSavedFields();

        //- Single phase flux
        void update
        (
            const surfaceScalarField& rhoOwn,
            const surfaceScalarField& rhoNei,
            const volVectorField& U,
            const volScalarField& e,
            const volScalarField& p,
            const volScalarField& c,
            surfaceScalarField& phi,
            surfaceScalarField& rhoPhi,
            surfaceVectorField& rhoUPhi,
            surfaceScalarField& rhoEPhi
        );

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

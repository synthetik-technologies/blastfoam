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
    Foam::shallowFluxScheme

Description
    Base class for flux schemes to interpolate fields and loop over faces
    and boundaries for a single shared velocity and energy

SourceFiles
    shallowFluxScheme.C
    shallowFluxSchemeNew.C

\*---------------------------------------------------------------------------*/

#ifndef shallowFluxScheme_H
#define shallowFluxScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluxSchemeBase.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class shallowFluxScheme Declaration
\*---------------------------------------------------------------------------*/

class shallowFluxScheme
:
    public fluxSchemeBase
{
protected:

    //- Dictionary for shallowFluxScheme
    const dictionary& dict_;

    //- Volumetric flux
    surfaceScalarField& phi_;

    //- Height flux
    surfaceScalarField& hPhi_;

    //- Hight momentum flux
    surfaceVectorField& hUPhi_;

    //- Gravitational acceleration
    const vector& g_;

    //- Magnitude of graviational acceleration
    const scalar magg_;

    //- Saved h0Own/h0Nei fields
    tmp<surfaceScalarField> th0Own_;
    tmp<surfaceScalarField> th0Nei_;

    // Protected Functions

        //- Return wavespeed
        inline scalar wavespeed(const scalar h) const
        {
            return sqrt(magg_*h);
        }

        //- Calcualte fluxes
        virtual void calculateFluxes
        (
            const scalar& hOwn, const scalar& hNei,
            const vector& UOwn, const vector& UNei,
            const vector& Sf,
            scalar& phi,
            scalar& hPhi,
            vector& hUPhi,
            const label facei, const label patchi = -1
        ) = 0;

        //- Calcualte flux given own/nei values
        virtual scalar calculateFlux
        (
            const scalar& fOwn, const scalar& fNei,
            const scalar& phi,
            const label facei, const label patchi = -1
        ) const
        {
            NotImplemented;
            return fOwn;
        }

        //- Interpolate field
        virtual scalar interpolate
        (
            const scalar& fOwn, const scalar& fNei,
            const label facei, const label patchi = -1
        ) const
        {
            NotImplemented;
            return fOwn;
        }

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
    TypeName("shallowFluxScheme");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            shallowFluxScheme,
            dictionary,
            (
                surfaceScalarField& phi,
                surfaceScalarField& hPhi,
                surfaceVectorField& hUPhi,
                const dimensionedVector& g
            ),
            (phi, hPhi, hUPhi, g)
        );

    // Constructor
    shallowFluxScheme
    (
        surfaceScalarField& phi,
        surfaceScalarField& hPhi,
        surfaceVectorField& hUPhi,
        const dimensionedVector& g
    );


    //- Destructor
    virtual ~shallowFluxScheme();


    // Selectors

        static autoPtr<shallowFluxScheme> New
        (
            surfaceScalarField& phi,
            surfaceScalarField& hPhi,
            surfaceVectorField& hUPhi,
            const dimensionedVector& g
        );



    // Member Functions

        //- Clear savedFields
        virtual void clear();

        //- Allocate saved fields
        virtual void createSavedFields();

        //- Single phase flux
        void update
        (
            const volScalarField& h,
            const volScalarField& h0,
            const volVectorField& U,
            const volVectorField& hU
        );

        //- Return the elevation source
        virtual tmp<volVectorField> ghGradH0
        (
            const dimensionedVector& g,
            const volScalarField& h,
            const volScalarField& h0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

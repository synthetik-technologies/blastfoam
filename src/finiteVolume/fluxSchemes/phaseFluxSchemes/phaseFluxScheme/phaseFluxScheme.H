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
    Foam::phaseFluxScheme

Description
    Base class for flux schemes to interpolate fields and loop over faces
    and boundaries for a phase with a unique velocity and energy

SourceFiles
    phaseFluxScheme.C
    phaseFluxSchemeNew.C

\*---------------------------------------------------------------------------*/

#ifndef phaseFluxScheme_H
#define phaseFluxScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluxSchemeBase.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseFluxScheme Declaration
\*---------------------------------------------------------------------------*/

class phaseFluxScheme
:
    public fluxSchemeBase
{
protected:

    //- Phase name
    const word phaseName_;

    //- Dictionary for phase fluxScheme
    const dictionary& dict_;

    //- Saved interpolated alpha field
    tmp<surfaceScalarField> alphaf_;

    //- Saved interpolated U field
    tmp<surfaceVectorField> Uf_;

    //- Saved interpolated p field
    tmp<surfaceScalarField> pf_;

    // Protected Functions

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
            scalar& alphaPhi,
            scalar& alphaRhoPhi,
            vector& alphaRhoUPhi,
            scalar& alphaRhoEPhi,
            const label facei, const label patchi = -1
        ) = 0;


        //- Update fields before calculating fluxes
        virtual void preUpdate(const volScalarField& p)
        {}

        //- Correct fluxes
        virtual void postUpdate()
        {}

        //- Return the name of a give field
        word fieldName(const word& name) const
        {
            return IOobject::groupName(type() + ":" + name, phaseName_);
        }


public:

    //- Runtime type information
    TypeName("phaseFluxScheme");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseFluxScheme,
            dictionary,
            (
                const surfaceScalarField& phi
            ),
            (phi)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseFluxScheme,
            solid,
            (
                const surfaceScalarField& phi
            ),
            (phi)
        );

    // Constructor
    phaseFluxScheme(const surfaceScalarField& phi);


    //- Destructor
    virtual ~phaseFluxScheme();


    // Selectors

        //- General selector
        static autoPtr<phaseFluxScheme> New
        (
            const surfaceScalarField& phi
        );

        //- Solid specific selection
        static autoPtr<phaseFluxScheme> NewSolid
        (
            const surfaceScalarField& phi
        );


    // Member Functions

        //- Clear savedFields
        virtual void clear();

        //- Allocate saved fields
        virtual void createSavedFields();

        //- Update
        void update
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& e,
            const volScalarField& p,
            const volScalarField& c,
            surfaceScalarField& phi,
            surfaceScalarField& alphaPhi,
            surfaceScalarField& alphaRhoPhi,
            surfaceVectorField& alphaRhoUPhi,
            surfaceScalarField& alphaRhoEPhi
        );

        //- Update
        void update
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& e,
            const volScalarField& p,
            const volScalarField& c,
            surfaceScalarField& phi,
            surfaceScalarField& alphaRhoPhi,
            surfaceVectorField& alphaRhoUPhi,
            surfaceScalarField& alphaRhoEPhi
        );

        //- Return interpolated U field
        tmp<surfaceVectorField> Uf() const;

        //- Return interpolated p field
        tmp<surfaceScalarField> pf() const;

        //- Return interpolated volume fraction field
        tmp<surfaceScalarField> alphaf() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

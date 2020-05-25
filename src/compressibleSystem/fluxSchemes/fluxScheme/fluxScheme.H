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
    and boundaries

SourceFiles
    fluxScheme.C
    newFluxScheme.C

\*---------------------------------------------------------------------------*/

#ifndef fluxScheme_H
#define fluxScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "surfaceFields.H"
#include "regIOobject.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "fvc.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluxScheme Declaration
\*---------------------------------------------------------------------------*/

class fluxScheme
:
    public regIOobject
{
protected:

    //- Const reference to mesh
    const fvMesh& mesh_;

    //- Owner and neighbour surface fields
    tmp<surfaceScalarField> own_;
    tmp<surfaceScalarField> nei_;

    //- Saved interpolated U field
    tmp<surfaceVectorField> Uf_;

    //- Saved interpolated rho fields
    tmp<surfaceScalarField> rhoOwn_;
    tmp<surfaceScalarField> rhoNei_;


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
            const label facei, const label patchi = -1
        ) const = 0;

        //- Interpolate field
        virtual scalar interpolate
        (
            const scalar& fOwn, const scalar& fNei,
            const label facei, const label patchi = -1
        ) const = 0;

        //- Update fields before calculating fluxes
        virtual void preUpdate(const volScalarField& p)
        {}

        //- Correct fluxes
        virtual void postUpdate()
        {}

        //- Return velocity w.r.t. face
        scalar meshPhi(const label facei, const label patchi) const
        {
            if (!mesh_.moving())
            {
                return 0.0;
            }
            if (patchi != -1)
            {
                return mesh_.phi().boundaryField()[patchi][facei];
            }
            else
            {
                return mesh_.phi()[facei];
            }
        }


        //- Set Uf at facei
        template<class Type>
        Type save
        (
            const label facei,
            const label patchi,
            const Type& x,
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& xf
        )
        {
            if (xf.valid())
            {
                if (patchi != -1)
                {
                    xf.ref().boundaryFieldRef()[patchi][facei] = x;
                }
                else
                {
                    xf.ref()[facei] = x;
                }
            }
            return x;
        }

        //- Set Uf at facei
        template<class Type>
        Type getValue
        (
            const label facei,
            const label patchi,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& xf
        ) const
        {
            if (patchi != -1)
            {
                return xf.boundaryField()[patchi][facei];
            }
            else
            {
                return xf[facei];
            }
        }

        inline word scheme(const word& name) const
        {
            word schemeName = "reconstruct(" + name + ")";
            if
            (
                !mesh_.schemesDict().subDict("interpolationSchemes").found(schemeName)
            )
            {
                WarningInFunction
                    << "Riemann schemes are used, but no limiter is " << nl
                    << "specified for " << name << "." << nl
                    << "This may result in unstable solutions." << endl;
            }
            return schemeName;
        }


public:

    //- Runtime type information
    TypeName("fluxScheme");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluxScheme,
            dictionary,
            (const fvMesh& mesh),
            (mesh)
        );

    // Constructor
    fluxScheme(const fvMesh& mesh);


    //- Destructor
    virtual ~fluxScheme();


    // Selectors

        static autoPtr<fluxScheme> New(const fvMesh& e);


    // Member Functions

        //- Clear savedFields
        virtual void clear();

        //- Allocate saved fields
        virtual void createSavedFields();

        //- Flux for three scalar fields
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& f,
            const word& fName
        ) const;

        //- Update
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

        //- Update
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

        //- Update
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

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fluxSchemeTmp.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

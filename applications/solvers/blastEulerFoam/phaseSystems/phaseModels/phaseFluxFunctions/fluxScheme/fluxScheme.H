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
        ) = 0;

        //- Interpolate field
        virtual scalar interpolate
        (
            const scalar& fOwn, const scalar& fNei,
            const bool rho,
            const label facei, const label patchi = -1
        ) const = 0;

        //- Update fields before calculating fluxes
        virtual void preUpdate(const volScalarField& p)
        {}

        //- Correct fluxes
        virtual void postUpdate()
        {}


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

        inline word scheme(const word& name, const word& phase) const
        {
            word schemeName =
                "reconstruct(" + IOobject::groupName(name, phase) + ")";
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
            (const fvMesh& mesh, const word& name),
            (mesh, name)
        );

    // Constructor
    fluxScheme(const fvMesh& mesh, const word& name);


    //- Destructor
    virtual ~fluxScheme();


    // Selectors

        static autoPtr<fluxScheme> New(const fvMesh& mesh, const word& name);


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

        //- Update
        void update
        (
            const UPtrList<volScalarField>& alphas,
            const UPtrList<volScalarField>& rhos,
            const volVectorField& U,
            const volScalarField& e,
            const volScalarField& p,
            const volScalarField& c,
            surfaceScalarField& phi,
            surfaceScalarField& alphaPhi,
            surfaceScalarField& alphaRhoPhi,
            PtrList<surfaceScalarField>& alphaPhis,
            PtrList<surfaceScalarField>& alphaRhoPhis,
            surfaceVectorField& alphaRhoUPhi,
            surfaceScalarField& alphaRhoEPhi
        );

        //- Return interpolated U field
        tmp<surfaceVectorField> Uf() const;

        //- Return interpolated p field
        tmp<surfaceScalarField> pf() const;

        //- Return interpolated volume fraction field
        tmp<surfaceScalarField> alphaf() const;

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

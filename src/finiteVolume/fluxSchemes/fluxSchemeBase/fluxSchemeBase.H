/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technologies
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
    Foam::fluxSchemeBase

Description
    Base class for flux schemes

SourceFiles
    fluxSchemeBase.C

\*---------------------------------------------------------------------------*/

#ifndef fluxSchemeBase_H
#define fluxSchemeBase_H

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
                           Class fluxSchemeBase Declaration
\*---------------------------------------------------------------------------*/

class fluxSchemeBase
:
    public regIOobject
{
public:

    static bool needEnergyFlux;

    //- Return the flux scheme
    static const fluxSchemeBase& findFluxScheme
    (
        const surfaceScalarField& phi
    );

    //- Does the flux scheme exist
    static bool foundFluxScheme
    (
        const surfaceScalarField& phi
    );

protected:
    // Protected data

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Volumetric flux
        const surfaceScalarField& phi_;

        //- Phases belonging to this flux scheme
        wordHashSet phases_;

    // Protected member functions

        //- Return velocity w.r.t. face
        scalar meshPhi(const label facei, const label patchi) const
        {
            if (!mesh_.moving())
            {
                return 0.0;
            }
            return
                patchi == -1
              ? mesh_.phi()[facei]
              : mesh_.phi().boundaryField()[patchi][facei];
        }


        //- Set a surface field at facei
        template<class Type>
        inline Type save
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
        inline Type getValue
        (
            const label facei,
            const label patchi,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& xf
        ) const
        {
            return
                patchi == -1
              ? xf[facei]
              : xf.boundaryField()[patchi][facei];
        }

        //- Set Uf at facei
        template<class Type>
        inline Type getValue
        (
            const label facei,
            const label patchi,
            const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& xf
        ) const
        {
            return
                patchi == -1
              ? xf()[facei]
              : xf().boundaryField()[patchi][facei];
        }


        //- Interpolate component of a field
        virtual scalar interpolate
        (
            const scalar& fOwn, const scalar& fNei,
            const label facei, const label patchi = -1
        ) const = 0;

        //- Compute Riemann flux of a given variable
        //  assumed to be a slave variable
        virtual scalar calculateFlux
        (
            const scalar& fOwn, const scalar& fNei,
            const scalar& phi,
            const label facei, const label patchi = -1
        ) const = 0;


public:

    //- Runtime type information
    TypeName("fluxScheme");

    // Constructor
    fluxSchemeBase(const surfaceScalarField& phi);


    //- Destructor
    virtual ~fluxSchemeBase();


    // Member Functions

        //- Return the volumetric flux
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Return the slave phases
        const wordHashSet& phases() const
        {
            return phases_;
        }

        wordHashSet& phases()
        {
            return phases_;
        }

        // Interpolation functions

            //- Return the upwind flux
            tmp<surfaceScalarField> upwindFlux() const;

            //- Templated interpolation
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
            (
                const GeometricField<Type, fvPatchField, volMesh>& f
            ) const;

            //- Templated interpolation
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& fOwn,
                const GeometricField<Type, fvsPatchField, surfaceMesh>& fNei,
                const word& name
            ) const;

            //- Templated interpolation
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
            (
                const GeometricField<Type, fvPatchField, volMesh>& f,
                const word& fName
            ) const;

            //- Templated interpolation for a phase, i.e. values at the interface
            //  are use upwinding based on interface direction
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> phaseInterpolate
            (
                const volScalarField& alpha,
                const GeometricField<Type, fvPatchField, volMesh>& f,
                const scalar rAlpha
            ) const;

            //- Templated interpolation for a phase
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> phaseInterpolate
            (
                const volScalarField& alpha,
                const GeometricField<Type, fvsPatchField, surfaceMesh>& fOwn,
                const GeometricField<Type, fvsPatchField, surfaceMesh>& fNei,
                const word& name,
                const scalar rAlpha
            ) const;

            //- Templated interpolation for a phase
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> phaseInterpolate
            (
                const volScalarField& alpha,
                const GeometricField<Type, fvPatchField, volMesh>& f,
                const word& fName,
                const scalar rAlpha
            ) const;


        // Fluxes

            //- Own/nei interpolation of f
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
            (
                const GeometricField<Type, fvPatchField, volMesh>& f,
                const surfaceScalarField& phi,
                const bool overwrite = true
            ) const;

            //- Own/nei interpolation of f and alphaRho
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
            (
                const GeometricField<Type, fvPatchField, volMesh>& f,
                const volScalarField& alphaRho,
                const surfaceScalarField& phi,
                const bool overwrite = true
            ) const;

            //- Own/nei interpolation of f
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
            (
                const GeometricField<Type, fvPatchField, volMesh>& f,
                const surfaceScalarField& alphaRhoOwn,
                const surfaceScalarField& alphaRhoNei,
                const surfaceScalarField& phi
            ) const;

            //- Own/nei interpolation of f
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& fOwn,
                const GeometricField<Type, fvsPatchField, surfaceMesh>& fNei,
                const surfaceScalarField& phi
            ) const;


            //- Own/nei phase interpolation of f
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> phaseFlux
            (
                const volScalarField& alpha,
                const GeometricField<Type, fvPatchField, volMesh>& f,
                const surfaceScalarField& phi,
                const scalar rAlpha,
                const bool overwrite = true
            ) const;

            //- Own/nei interpolation of f and alphaRho
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> phaseFlux
            (
                const volScalarField& alpha,
                const GeometricField<Type, fvPatchField, volMesh>& f,
                const volScalarField& alphaRho,
                const surfaceScalarField& phi,
                const scalar rAlpha,
                const bool overwrite = true
            ) const;

            //- Own/nei interpolation of f
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> phaseFlux
            (
                const volScalarField& alpha,
                const GeometricField<Type, fvPatchField, volMesh>& f,
                const surfaceScalarField& alphaRhoOwn,
                const surfaceScalarField& alphaRhoNei,
                const surfaceScalarField& phi,
                const scalar rAlpha
            ) const;

            //- Own/nei interpolation of f
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> phaseFlux
            (
                const volScalarField& alpha,
                const GeometricField<Type, fvsPatchField, surfaceMesh>& fOwn,
                const GeometricField<Type, fvsPatchField, surfaceMesh>& fNei,
                const surfaceScalarField& phi,
                const scalar rAlpha
            ) const;


        //- Correct phase fields
        template<class Type>
        static void correctPhaseFields
        (
            const volScalarField& alpha,
            GeometricField<Type, fvsPatchField, surfaceMesh>& fOwn,
            GeometricField<Type, fvsPatchField, surfaceMesh>& fNei,
            const scalar rAlpha
        );

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fluxSchemeBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

protected:
    // Protected data

        //- Const reference to mesh
        const fvMesh& mesh_;


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
            const bool rho,
            const label facei, const label patchi = -1
        ) const = 0;


public:

    //- Runtime type information
    TypeName("fluxScheme");

    // Constructor
    fluxSchemeBase(const fvMesh& mesh, const word& phaseName = word::null);


    //- Destructor
    virtual ~fluxSchemeBase();


    // Member Functions

        //- Interpolation functions
        virtual tmp<surfaceScalarField> interpolate
        (
            const volScalarField& f,
            const word& fName
        ) const;
        virtual tmp<surfaceVectorField> interpolate
        (
            const volVectorField& f,
            const word& fName
        ) const;
        virtual tmp<surfaceSymmTensorField> interpolate
        (
            const volSymmTensorField& f,
            const word& fName
        ) const;
        virtual tmp<surfaceSphericalTensorField> interpolate
        (
            const volSphericalTensorField& f,
            const word& fName
        ) const;
        virtual tmp<surfaceTensorField> interpolate
        (
            const volTensorField& f,
            const word& fName
        ) const;

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
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolateField
        (
            const GeometricField<Type, fvPatchField, volMesh>& f,
            const word& fName
        ) const;


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

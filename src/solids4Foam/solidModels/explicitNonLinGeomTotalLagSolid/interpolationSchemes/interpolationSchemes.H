/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::interpolationSchemes

Description
    Define interpolation schemes.

SourceFiles
    interpolationSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef interpolationSchemes_H
#define interpolationSchemes_H

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "pointConstraints.H"
#include "volPointInterpolation.H"
#include "fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class interpolationSchemes
{
    // Private data

        //- Mesh
        const fvMesh& mesh_;

        //- Cell owners
        const labelUList& own_;

        //- Material cell centre coordinates
        const volVectorField& C_;

        //- Material face centre coordinates
        const surfaceVectorField& Cf_;

        //- Material nodal coordinates
        const pointField& points_;


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        interpolationSchemes(const interpolationSchemes&);

        //- Disallow default bitwise assignment
        void operator=(const interpolationSchemes&);

public:

    //- Runtime type information
    TypeName("interpolationSchemes");

    // Constructors

        //- Construct from mesh
            explicit interpolationSchemes(const fvMesh&);

    //- Destructor
    virtual ~interpolationSchemes();


    // Member functions

        // Edit

            //- Surface to volume interpolation
            tmp<volVectorField> surfaceToVol
            (
                const GeometricField<vector, fvsPatchField, surfaceMesh>&,
                const GeometricField<vector, pointPatchField, pointMesh>&
            ) const;

            //- Helper: push master point data to collocated points
            //  (Borrowed from volPointInterpolation class in OpenFOAM)
            template<class Type>
            void pushUntransformedData(List<Type>&) const;

            //- Add separated contributions
            //  (Borrowed from volPointInterpolation class in OpenFOAM)
            template<class Type>
            void addSeparated
            (
                GeometricField<Type, pointPatchField, pointMesh>&
            ) const;

            // Volume to point interpolation for vectorField
            void volToPoint
            (
                const GeometricField<vector, fvPatchField, volMesh>& U,
                const GeometricField<tensor, fvPatchField, volMesh>& gradU,
                GeometricField<vector, pointPatchField, pointMesh>& pU
            ) const;

            //- Point to surface interpolation
            tmp<surfaceVectorField> pointToSurface
            (
                const GeometricField<vector, pointPatchField, pointMesh>& U
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

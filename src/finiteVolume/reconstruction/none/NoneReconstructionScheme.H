/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020
     \\/     M anipulation  | Synthetik Applied Technology
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::NoneReconstructionScheme

Description
    Replaces MUSCL reconstruction with standard interpolation

SourceFiles
    NoneReconstructionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef NoneReconstructionScheme_H
#define NoneReconstructionScheme_H

#include "ReconstructionScheme.H"
#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class NoneReconstructionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class NoneReconstructionScheme
:
    public ReconstructionScheme<Type>
{
    //- Surface interpolation scheme
    tmp<surfaceInterpolationScheme<Type>> interp_;


protected:
// Protected data

    //- Calculate the limiter
    virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    calcLimiter(const scalar& dir) const
    {
        NotImplemented;
        return tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>();
    }


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from mesh interpolation scheme name (provided in New)
        NoneReconstructionScheme
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            Istream& is
        );

    //- Destructor
    virtual ~NoneReconstructionScheme();

    // Member Functions

        //- Return the owner interpolated field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolateOwn() const;

        //- Return the neighbor interpolated field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolateNei() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoneReconstructionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

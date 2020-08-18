/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
18-08-2020 Jeff Heylmun:    | Use standard interpolation
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
    Foam::noneMUSCLReconstructionScheme

Description
    Replaces MUSCL reconstruction with standard interpolation

SourceFiles
    noneMUSCLReconstructionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef noneMUSCLReconstructionScheme_H
#define noneMUSCLReconstructionScheme_H

#include "MUSCLReconstructionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class noneMUSCLReconstructionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class noneMUSCLReconstructionScheme
:
    public MUSCLReconstructionScheme<Type>
{
protected:
// Protected data

    //- Name of interpolation scheme
    word name_;

    //- Owner interpolation field
    const surfaceScalarField& own_;

    //- Neighbour interpolation field
    const surfaceScalarField& nei_;

    //- Return reference to field if found in db
    const surfaceScalarField& lookupOrConstruct
    (
        const word& fieldName,
        const scalar& value
    ) const;

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
        noneMUSCLReconstructionScheme
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            Istream& is
        );

    //- Destructor
    virtual ~noneMUSCLReconstructionScheme();

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
    #include "noneMUSCLReconstructionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
18-08-2020 Jeff Heylmun:    | Quadratic interpolation with MUSCL reconstruction
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
    Foam::quadraticMUSCLReconstructionScheme

Description
    Quadratic MUSCL reconstruction using gradient and Hessian of field
    components. Interpolated values are limited to the min/max of the own/nei
    cell values.


SourceFiles
    quadraticMUSCLReconstructionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticMUSCLReconstructionScheme_H
#define quadraticMUSCLReconstructionScheme_H

#include "MUSCLReconstructionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class KimScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class quadraticMUSCLReconstructionScheme
:
    public MUSCLReconstructionScheme<Type>
{
protected:
// Protected data

    //- Gradient of field
    PtrList<GeometricField<vector, fvPatchField, volMesh>> gradPhis_;

    //- Hessian of field
    PtrList<GeometricField<tensor, fvPatchField, volMesh>> hessPhis_;


public:

    //- Runtime type information
    TypeName("quadraticMUSCL");


    // Constructors

        //- Construct from mesh and faceFlux and limiter scheme
        quadraticMUSCLReconstructionScheme
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            Istream& is
        );

    //- Destructor
    virtual ~quadraticMUSCLReconstructionScheme();

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
    #include "quadraticMUSCLReconstructionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

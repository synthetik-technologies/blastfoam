/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
18-08-2020 Jeff Heylmun:    | Upwind MUSCL reconstruction
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
    Foam::upwindMUSCLReconstructionScheme

Description
    Owner and neighbor interpolated fields are set to the owner and neighbor
    cell values.

SourceFiles
    upwindMUSCLReconstructionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef upwindMUSCLReconstructionScheme_H
#define upwindMUSCLReconstructionScheme_H

#include "MUSCLReconstructionScheme.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class KimScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class upwindMUSCLReconstructionScheme
:
    public MUSCLReconstructionScheme<Type>
{

     //- Calculate the limiter
    virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    calcLimiter(const scalar& dir) const
    {
        NotImplemented;
        return tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>();
    }


public:

    //- Runtime type information
    TypeName("upwindMUSCL");


    // Constructors

        //- Construct from mesh and faceFlux and limiter scheme
        upwindMUSCLReconstructionScheme
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            Istream& is
        );

    //- Destructor
    virtual ~upwindMUSCLReconstructionScheme();


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
    #include "upwindMUSCLReconstructionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

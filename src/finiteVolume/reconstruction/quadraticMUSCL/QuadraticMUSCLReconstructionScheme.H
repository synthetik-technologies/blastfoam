/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2024
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
    Foam::QuadraticMUSCLReconstructionScheme

Description
    Quadratic MUSCL reconstruction using gradient and Hessian of field
    components. Interpolated values are limited to the min/max of the own/nei
    cell values.


SourceFiles
    QuadraticMUSCLReconstructionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef QuadraticMUSCLReconstructionScheme_H
#define QuadraticMUSCLReconstructionScheme_H

#include "ReconstructionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class KimScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class QuadraticMUSCLReconstructionScheme
:
    public ReconstructionScheme<Type>
{
protected:
// Protected data

    //- Gradient of field
    mutable PtrList<GeometricField<vector, fvPatchField, volMesh>> gradPhis_;

    //- Hessian of field
    mutable PtrList<GeometricField<tensor, fvPatchField, volMesh>> hessPhis_;

    //- Bound face values
    bool bound_;

    //- Extrapolate boundaries
    bool extrapolate_;

    inline tmp<volVectorField> gradient
    (
        const volScalarField&,
        const label diri
    ) const
    {
        return gradPhis_[diri];
    }

    //- Read optional inputs
    void read(Istream& is)
    {
        if (is.good())
        {
            bound_ = readBool(is);
        }
        if (is.good())
        {
            extrapolate_ = readBool(is);
        }
    }

    //- Construct gradients
    void constructGradPhis() const;


public:

    //- Runtime type information
    TypeName("quadraticMUSCL");


    // Constructors

        //- Construct from mesh and faceFlux and limiter scheme
        QuadraticMUSCLReconstructionScheme
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            Istream& is,
            const bool overwrite
        );

    //- Destructor
    virtual ~QuadraticMUSCLReconstructionScheme();

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
    #include "QuadraticMUSCLReconstructionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

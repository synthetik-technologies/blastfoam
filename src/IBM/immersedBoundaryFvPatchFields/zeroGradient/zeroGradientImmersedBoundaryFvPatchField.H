/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
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
    Foam::zeroGradientImmersedBoundaryFvPatchField

Description


SourceFiles
    zeroGradientImmersedBoundaryFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef zeroGradientImmersedBoundaryFvPatchField_H
#define zeroGradientImmersedBoundaryFvPatchField_H

#include "volFields.H"
#include "immersedBoundaryFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class zeroGradientImmersedBoundaryFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class zeroGradientImmersedBoundaryFvPatchField
:
    public immersedBoundaryFvPatchField<Type>
{
    //- Number of smoothing iterations
    label nSmooth_;

public:

    //- Runtime type information
    TypeName("zeroGradient");

    // Constructor
    zeroGradientImmersedBoundaryFvPatchField
    (
        GeometricField<Type, fvPatchField, volMesh>&,
        const dictionary&,
        const immersedBoundaryObject&
    );

    // Destructor
    virtual ~zeroGradientImmersedBoundaryFvPatchField()
    {}


    // Member Functions


        // Evaluation functions

            //- Add forcing to the given field
            virtual void addForcing
            (
                Field<Type>& f,
                const Field<scalar>& alphaRho,
                const Field<Type>& old,
                const Field<Type>& RHS,
                const scalar& dt
            ) const;

            //- Update size and values
            virtual void updateCoeffs() const;

            //- Add forcing to the given field
            virtual void setValues();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zeroGradientImmersedBoundaryFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Description

Usage

See also
    Foam::couple

SourceFiles
    turbulentTemperatureRadCoupledBlastMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef solidTractionCoupledImmersedFvPatchVectorField_H
#define solidTractionCoupledImmersedFvPatchVectorField_H

#include "fixedValueFvPatchField.H"
#include "immersedMeshMapper.H"
#include "immersedBoundaryObjectListSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class turbulentTemperatureRadCoupledBlastMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class coupledImmersedFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private Data

        //- Name of object
        const word name_;

        //- Mapper
        const immersedMeshMapper& mapper_;

        //- Immersed object
        const immersedBoundaryObject& ibm_;

public:

    //- Runtime type information
    TypeName("coupledImmersed");


    // Constructors

        //- Construct from patch and internal field
        coupledImmersedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        coupledImmersedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  coupledImmersedFvPatchField onto a
        //  new patch
        coupledImmersedFvPatchField
        (
            const coupledImmersedFvPatchField&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        coupledImmersedFvPatchField
        (
            const coupledImmersedFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new coupledImmersedFvPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "coupledImmersedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

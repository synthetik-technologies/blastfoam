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
    Foam::coupledPointMotionFvPatchField

Description
    Foam::coupledMotionPointPatchField

SourceFiles
    coupledMotionPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledMotionPointPatchField_H
#define coupledMotionPointPatchField_H

#include "Random.H"
#include "fixedValuePointPatchFields.H"
#include "mappedPointPatchSelector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class coupledMotionPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class coupledMotionPointPatchField
:
    public fixedValuePointPatchField<Type>
{
    // Private Data

        mappedPointPatchSelector mpp_;

        //- Name of neighbour displacement field
        const word DnbrName_;

        //- The component to solve for
        direction cmpt_;

        //- Is the field a velocity
        const bool velocity_;

        direction cmpt(const word& cmptName) const
        {
            if (cmptName == "x")
            {
                return vector::X;
            }
            else if (cmptName == "y")
            {
                return vector::Y;
            }
            else if (cmptName == "z")
            {
                return vector::Z;
            }
            else
            {
                FatalErrorInFunction
                    << "Given component name " << cmptName << " should be x, y or z"
                    << exit(FatalError);

                return 0;
            }
        }

        tmp<Field<Type>> getCmpt(const vectorField& f) const;


public:

    //- Runtime type information
    TypeName("coupledPointMotion");


    // Constructors

        //- Construct from patch and internal field
        coupledMotionPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        coupledMotionPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given coupledMotionPointPatchField
        //  onto a new patch
        coupledMotionPointPatchField
        (
            const coupledMotionPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Copy constructor
        coupledMotionPointPatchField
        (
            const coupledMotionPointPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        coupledMotionPointPatchField
        (
            const coupledMotionPointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new coupledMotionPointPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "coupledMotionPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

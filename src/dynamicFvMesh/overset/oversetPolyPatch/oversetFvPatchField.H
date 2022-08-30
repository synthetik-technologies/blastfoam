/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::oversetFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    Boundary condition for use on overset patches. To be run in combination
    with special dynamicFvMesh type that inserts interpolation into the matrix.

SourceFiles
    oversetFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef oversetFvPatchField_H
#define oversetFvPatchField_H

#include "oversetFvPatch.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class oversetFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class oversetFvPatchField
:
    public zeroGradientFvPatchField<Type>
{
protected:

    // Protected data

        //- Local reference cast into the overset patch
        const oversetFvPatch& oversetPatch_;

        //- Master patch ID
        mutable label masterPatchID_;


public:

    //- Runtime type information
    TypeName(oversetFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        oversetFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        oversetFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given oversetFvPatchField onto a new patch
        oversetFvPatchField
        (
            const oversetFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        oversetFvPatchField(const oversetFvPatchField<Type>&) = delete;

        //- Construct as copy setting internal field reference
        oversetFvPatchField
        (
            const oversetFvPatchField<Type>&,
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
                new oversetFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "oversetFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

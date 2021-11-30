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
    Foam::wedgeLinearMomentumFvPatchVectorField

Description
    Specifies contact linear momentum for a wedge boundary patch
    according to the Riemann solver. External boundary traction can
    be specified together with a linear ramp function.

SeeAlso
    Foam::wedgeFvPatchField

SourceFiles
    wedgeLinearMomentumFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef wedgeLinearMomentumFvPatchVectorField_H
#define wedgeLinearMomentumFvPatchVectorField_H

#include "wedgeFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class wedgeLinearMomentumFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class wedgeLinearMomentumFvPatchVectorField
:
    public wedgeFvPatchVectorField
{
public:

    // Runtime type information
    TypeName("wedgeLinearMomentum");


    // Constructors

        //- Construct from patch and internal field
        wedgeLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wedgeLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  wedgeLinearMomentumFvPatchVectorField onto a new patch
        wedgeLinearMomentumFvPatchVectorField
        (
            const wedgeLinearMomentumFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wedgeLinearMomentumFvPatchVectorField
        (
            const wedgeLinearMomentumFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        wedgeLinearMomentumFvPatchVectorField
        (
            const wedgeLinearMomentumFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new wedgeLinearMomentumFvPatchVectorField(*this, iF)
            );
        }

    // Member functions

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        // Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

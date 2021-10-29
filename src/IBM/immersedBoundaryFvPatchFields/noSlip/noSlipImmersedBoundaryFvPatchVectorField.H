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
    Foam::noSlipImmersedBoundaryFvPatchField

Description
    This boundary condition is not designed to be evaluated; it is assmued
    that the value is assigned via field assignment, and not via a call to
    e.g. \c updateCoeffs or \c evaluate.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            noSlipImmersedBoundary;
    }
    \endverbatim

SourceFiles
    noSlipImmersedBoundaryFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef noSlipImmersedBoundaryFvPatchVectorField_H
#define noSlipImmersedBoundaryFvPatchVectorField_H

#include "immersedBoundaryFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class noSlipImmersedBoundaryFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class noSlipImmersedBoundaryFvPatchVectorField
:
    public immersedBoundaryVectorPatchField
{
public:

    //- Runtime type information
    TypeName("noSlip");


    // Constructors

        //- Construct from patch and internal field
        noSlipImmersedBoundaryFvPatchVectorField
        (
            volVectorField&,
            const dictionary&,
            const immersedBoundaryObject&
        );


    // Member Functions


        // Evaluation functions

            //- Update size and values
            virtual void updateCoeffs() const;

            //- Add forcing to the given field
            virtual void setValues();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
